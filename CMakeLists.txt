cmake_minimum_required(VERSION 3.15)
project(bpm)

include(FetchContent)

FetchContent_Declare(
  libsais
  GIT_REPOSITORY https://github.com/IlyaGrebnov/libsais.git
  GIT_TAG master
)

FetchContent_Declare(
  sdsl
  GIT_REPOSITORY https://github.com/simongog/sdsl-lite.git
  GIT_TAG master
)

FetchContent_MakeAvailable(libsais sdsl)

add_executable(bpm src/main.cpp)
target_link_libraries(bpm PRIVATE libsais sdsl)
