cmake_minimum_required(VERSION 3.10)

project(minimum_weight_and_diameter_spanning_tree VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(mvdst ${PROJECT_SOURCES}
  src/main.cpp
  src/mvdst.cpp
  src/io.cpp
  src/vertex_impl.cpp
  )

target_include_directories(mvdst PRIVATE src)

target_include_directories(mvdst PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/lib/boost/boost/
  src/
  )

set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -O3")

