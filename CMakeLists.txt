cmake_minimum_required(VERSION 3.15)
set(CMAKE_POLICY_VERSION_MINIMUM 3.15)
cmake_policy(SET CMP0148 NEW)

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
    GIT_REPOSITORY https://github.com/vgteam/sdsl-lite.git
    GIT_TAG master
)

FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG master
)

FetchContent_MakeAvailable(libsais sdsl doctest)

add_library(bpm_core
    src/motifs_search.cpp
    src/preprocessing.cpp
    src/data_import.cpp
)

target_include_directories(bpm_core PUBLIC
    src
    ${libsais_SOURCE_DIR}
    ${sdsl_SOURCE_DIR}/include
    ${sdsl_BINARY_DIR}/external/sdsl-lite/include
)
target_link_libraries(bpm_core PUBLIC libsais sdsl)

add_executable(bpm src/main.cpp)
target_link_libraries(bpm PRIVATE bpm_core)

add_executable(run_tests tests/test_main.cpp)
target_link_libraries(run_tests PRIVATE bpm_core doctest::doctest)
