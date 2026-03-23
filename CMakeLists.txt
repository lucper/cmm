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
    GIT_TAG        master
)
FetchContent_MakeAvailable(libsais)

FetchContent_Declare(
    gtl
    GIT_REPOSITORY https://github.com/greg7mdp/gtl.git
    GIT_TAG        v1.2.0
)
FetchContent_MakeAvailable(gtl)

add_library(bpm_core
    src/motifs_search.cpp
    src/rank_table.cpp
    src/data_import.cpp
)
target_link_libraries(bpm_core PUBLIC libsais gtl)
target_compile_options(bpm_core PRIVATE -g -fno-omit-frame-pointer) # profiling

add_executable(bpm src/main.cpp)
target_link_libraries(bpm PRIVATE bpm_core)
target_compile_options(bpm PRIVATE -fno-omit-frame-pointer) # profiling

add_executable(run_tests tests/test_main.cpp)
target_link_libraries(run_tests PRIVATE bpm_core)
