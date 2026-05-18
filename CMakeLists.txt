cmake_minimum_required(VERSION 3.15)
set(CMAKE_POLICY_VERSION_MINIMUM 3.15)
cmake_policy(SET CMP0148 NEW)

project(cmm)

if(DEFINED USE_64)
    add_compile_definitions(USE_64)
elseif(DEFINED USE_32)
    add_compile_definitions(USE_32)
else()
    message(FATAL_ERROR "Set either -DUSE_32 or -DUSE_64.")
endif()

include(FetchContent)

FetchContent_Declare(
    libsais
    GIT_REPOSITORY https://github.com/IlyaGrebnov/libsais.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(libsais)

FetchContent_Declare(
    cxxopts
    GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(cxxopts)

FetchContent_Declare(
    indicators
    GIT_REPOSITORY https://github.com/p-ranav/indicators.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(indicators)

find_package(OpenMP REQUIRED)

add_library(cmm_core
    src/motifs_search.cpp
    src/rank_table.cpp
    src/data_import.cpp
)
target_link_libraries(cmm_core PUBLIC indicators cxxopts libsais OpenMP::OpenMP_CXX)
target_compile_options(cmm_core PRIVATE
    -Wall -Wextra
    -fno-omit-frame-pointer
    -O3 -march=native
    $<$<CONFIG:Release>:-flto>
    $<$<CONFIG:Debug>:-g>
)

add_executable(cmm src/main.cpp)
target_link_libraries(cmm PRIVATE cmm_core OpenMP::OpenMP_CXX)
target_compile_options(cmm PRIVATE
    $<$<CONFIG:Release>:-flto>
)

add_executable(cmm_perf experiments/performance/perf_main_algo.cpp)
target_include_directories(cmm_perf PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(cmm_perf PRIVATE cmm_core OpenMP::OpenMP_CXX)
target_compile_options(cmm_perf PRIVATE
    $<$<CONFIG:Release>:-flto>
)
