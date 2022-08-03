cmake_minimum_required(VERSION 3.12)
project(test-find-package-sais)

# You should set the CMAKE_PREFIX_PATH to the installed directory with libsais to find it,
# if it is not installed in your prefixes...
find_package(sais REQUIRED)

add_library(test-package main.cpp)
target_link_libraries(test-package sais::sais)
