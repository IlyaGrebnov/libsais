cmake_minimum_required(VERSION 3.12)
project(libsais VERSION 2.7.1 LANGUAGES C DESCRIPTION "libsais is a library for linear time suffix array, longest common prefix array and burrows wheeler transform construction based on induced sorting algorithm.")

option(BUILD_SHARED_LIBS "Build libsais as shared lib" OFF)

# For proper install in system...
include(GNUInstallDirs)

# The whole library with adjusted include_directories for `make install` users. They must use e.g. `#include <libsais/libsais.h>`.
add_library(sais src/libsais.c src/libsais.h src/libsais16.c src/libsais16.h src/libsais64.c src/libsais64.h)
target_include_directories(sais PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>" "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
set_target_properties(sais PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR})

# Installing target artifacts (libs), all targets must be inside an export set (in our case it's the sais-config)...
install(TARGETS sais EXPORT sais-config
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
# public headers..
install(FILES src/libsais.h src/libsais16.h src/libsais64.h DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/libsais")
# generate config with imported targets for find_package(sais) with generating the pretty namespaced targets like sais::sais for linking...
install(EXPORT sais-config NAMESPACE sais:: DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/sais")

# Helper for using the library without installation (via `make install`, e.g., include("thirdparty/libsais-2.7.1/build/sais-config.cmake"))
export(EXPORT sais-config NAMESPACE sais::)

# Helper for versioning the packages and libraries.. we assume that library is compatible with any minor versions...
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/sais-config-version.cmake" COMPATIBILITY SameMinorVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/sais-config-version.cmake" DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/sais")

# Install docs, LICENSE, etc...
install(FILES CHANGES LICENSE README.md VERSION DESTINATION "${CMAKE_INSTALL_DOCDIR}")

