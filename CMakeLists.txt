cmake_minimum_required(VERSION 3.16)
project(min_surf_netw_cpp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_compile_options(-Wall -Wextra)

add_executable(min-surf-netw
  src/main.cpp
  src/preset.cpp
  src/network.cpp
)
target_include_directories(min-surf-netw PRIVATE src)
