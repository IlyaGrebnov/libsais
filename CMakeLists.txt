cmake_minimum_required(VERSION 3.10)

project(libsais VERSION 2.9.1 LANGUAGES C DESCRIPTION "libsais is a library for linear time suffix array, longest common prefix array and burrows wheeler transform construction based on induced sorting algorithm.")

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(LIBSAIS_USE_OPENMP "Use OpenMP for parallelization" OFF)
option(LIBSAIS_BUILD_SHARED_LIB "Build libsais as a shared library" OFF)

if(LIBSAIS_BUILD_SHARED_LIB)
    set(LIBSAIS_LIBRARY_TYPE SHARED)
else()
    set(LIBSAIS_LIBRARY_TYPE STATIC)
endif()

add_library(libsais ${LIBSAIS_LIBRARY_TYPE})

target_sources(libsais PRIVATE
    include/libsais.h
    include/libsais16.h
    include/libsais16x64.h
    include/libsais64.h
    src/libsais.c
    src/libsais16.c
    src/libsais16x64.c
    src/libsais64.c
)

if(LIBSAIS_USE_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

if(LIBSAIS_USE_OPENMP AND OpenMP_C_FOUND)
    target_compile_definitions(libsais PUBLIC LIBSAIS_OPENMP)
    target_link_libraries(libsais PRIVATE OpenMP::OpenMP_C)
endif()

if(LIBSAIS_BUILD_SHARED_LIB)
    target_compile_definitions(libsais PUBLIC LIBSAIS_SHARED)
    target_compile_definitions(libsais PRIVATE LIBSAIS_EXPORTS)
endif()

target_include_directories(libsais PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
