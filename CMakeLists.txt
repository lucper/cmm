cmake_minimum_required(VERSION 3.15)
project(bpm)

if(DEFINED _USE_64)
    add_compile_definitions(_USE_64)
elseif(DEFINED _USE_32)
    add_compile_definitions(_USE_32)
else()
    message(FATAL_ERROR "Set either -D_USE_32 or -D_USE_64.")
endif()

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

add_executable(bpm src/main.cpp src/preprocessing.hpp src/preprocessing.cpp)
target_link_libraries(bpm PRIVATE libsais sdsl)
