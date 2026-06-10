// Minimal self-contained JSON parser/serializer (no external dependencies).
// Supports objects, arrays, strings, numbers, booleans and null, which is all
// the preset format requires. Not a general-purpose library: it aims to be
// small, correct for well-formed input, and dependency-free.
#pragma once

#include <cctype>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace msn {

class Json {
public:
    enum class Type { Null, Bool, Number, String, Array, Object };

    Type type = Type::Null;
    bool boolValue = false;
    double numberValue = 0.0;
    std::string stringValue;
    std::vector<Json> arrayValue;
    std::map<std::string, Json> objectValue;

    Json() = default;

    bool isObject() const { return type == Type::Object; }
    bool isArray() const { return type == Type::Array; }
    bool isNumber() const { return type == Type::Number; }

    bool has(const std::string& key) const {
        return type == Type::Object && objectValue.count(key) > 0;
    }
    const Json& at(const std::string& key) const {
        auto it = objectValue.find(key);
        if (it == objectValue.end()) throw std::runtime_error("JSON key not found: " + key);
        return it->second;
    }
    const Json& operator[](std::size_t i) const { return arrayValue[i]; }
    std::size_t size() const { return arrayValue.size(); }

    double asNumber() const {
        if (type != Type::Number) throw std::runtime_error("JSON value is not a number");
        return numberValue;
    }
    int asInt() const { return static_cast<int>(std::lround(asNumber())); }
    bool asBool() const {
        if (type == Type::Bool) return boolValue;
        if (type == Type::Number) return numberValue != 0.0;
        throw std::runtime_error("JSON value is not a bool");
    }
    const std::string& asString() const {
        if (type != Type::String) throw std::runtime_error("JSON value is not a string");
        return stringValue;
    }

    static Json parse(const std::string& text) {
        std::size_t pos = 0;
        Json j = parseValue(text, pos);
        skipWs(text, pos);
        if (pos != text.size()) throw std::runtime_error("Trailing characters in JSON");
        return j;
    }

private:
    static void skipWs(const std::string& s, std::size_t& p) {
        while (p < s.size() && (s[p] == ' ' || s[p] == '\t' || s[p] == '\n' || s[p] == '\r')) ++p;
        // Allow // line comments and /* */ block comments for human-edited presets.
        if (p + 1 < s.size() && s[p] == '/' && s[p + 1] == '/') {
            while (p < s.size() && s[p] != '\n') ++p;
            skipWs(s, p);
        } else if (p + 1 < s.size() && s[p] == '/' && s[p + 1] == '*') {
            p += 2;
            while (p + 1 < s.size() && !(s[p] == '*' && s[p + 1] == '/')) ++p;
            p += 2;
            skipWs(s, p);
        }
    }

    static Json parseValue(const std::string& s, std::size_t& p) {
        skipWs(s, p);
        if (p >= s.size()) throw std::runtime_error("Unexpected end of JSON");
        char c = s[p];
        if (c == '{') return parseObject(s, p);
        if (c == '[') return parseArray(s, p);
        if (c == '"') return parseString(s, p);
        if (c == 't' || c == 'f') return parseBool(s, p);
        if (c == 'n') { p += 4; Json j; j.type = Type::Null; return j; }
        return parseNumber(s, p);
    }

    static Json parseObject(const std::string& s, std::size_t& p) {
        Json j;
        j.type = Type::Object;
        ++p;  // {
        skipWs(s, p);
        if (p < s.size() && s[p] == '}') { ++p; return j; }
        while (true) {
            skipWs(s, p);
            Json key = parseString(s, p);
            skipWs(s, p);
            if (s[p] != ':') throw std::runtime_error("Expected ':' in JSON object");
            ++p;
            Json value = parseValue(s, p);
            j.objectValue[key.stringValue] = std::move(value);
            skipWs(s, p);
            if (p >= s.size()) throw std::runtime_error("Unterminated JSON object");
            if (s[p] == ',') { ++p; continue; }
            if (s[p] == '}') { ++p; break; }
            throw std::runtime_error("Expected ',' or '}' in JSON object");
        }
        return j;
    }

    static Json parseArray(const std::string& s, std::size_t& p) {
        Json j;
        j.type = Type::Array;
        ++p;  // [
        skipWs(s, p);
        if (p < s.size() && s[p] == ']') { ++p; return j; }
        while (true) {
            Json value = parseValue(s, p);
            j.arrayValue.push_back(std::move(value));
            skipWs(s, p);
            if (p >= s.size()) throw std::runtime_error("Unterminated JSON array");
            if (s[p] == ',') { ++p; continue; }
            if (s[p] == ']') { ++p; break; }
            throw std::runtime_error("Expected ',' or ']' in JSON array");
        }
        return j;
    }

    static Json parseString(const std::string& s, std::size_t& p) {
        Json j;
        j.type = Type::String;
        if (s[p] != '"') throw std::runtime_error("Expected '\"' in JSON");
        ++p;
        std::string out;
        while (p < s.size() && s[p] != '"') {
            char c = s[p++];
            if (c == '\\' && p < s.size()) {
                char e = s[p++];
                switch (e) {
                    case 'n': out.push_back('\n'); break;
                    case 't': out.push_back('\t'); break;
                    case 'r': out.push_back('\r'); break;
                    case '"': out.push_back('"'); break;
                    case '\\': out.push_back('\\'); break;
                    case '/': out.push_back('/'); break;
                    default: out.push_back(e); break;
                }
            } else {
                out.push_back(c);
            }
        }
        if (p >= s.size()) throw std::runtime_error("Unterminated JSON string");
        ++p;  // closing quote
        j.stringValue = std::move(out);
        return j;
    }

    static Json parseBool(const std::string& s, std::size_t& p) {
        Json j;
        j.type = Type::Bool;
        if (s.compare(p, 4, "true") == 0) { j.boolValue = true; p += 4; }
        else if (s.compare(p, 5, "false") == 0) { j.boolValue = false; p += 5; }
        else throw std::runtime_error("Invalid JSON literal");
        return j;
    }

    static Json parseNumber(const std::string& s, std::size_t& p) {
        std::size_t start = p;
        while (p < s.size() && (std::isdigit(static_cast<unsigned char>(s[p])) || s[p] == '+' ||
                                s[p] == '-' || s[p] == '.' || s[p] == 'e' || s[p] == 'E'))
            ++p;
        Json j;
        j.type = Type::Number;
        j.numberValue = std::stod(s.substr(start, p - start));
        return j;
    }
};

}  // namespace msn
