cmake_minimum_required(VERSION 3.16)
project(boris-pusher-box LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall -Wextra)

find_package(SFML 2.6 COMPONENTS graphics window system REQUIRED)

add_executable(boris-box
  src/main.cpp
  src/boris.cpp
  src/fields.cpp
)

target_link_libraries(boris-box PRIVATE sfml-graphics sfml-window sfml-system)
