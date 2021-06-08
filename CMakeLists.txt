# Copyright (c) 2017 Tyler Olsen
# Copyright (c) 2018-2019 Patrick Diehl
# Copyright (c) 2019 Nanmiao Wu
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.5)

project(BlazeIterative)

set(CMAKE_CXX_STANDARD 14)

#==========================================
# Optionally build tests (For devs)
#==========================================
option(BUILD_TESTS "Build BlazeIterative tests" OFF)

#==========================================
# Set up BlazeIterative target
#==========================================
add_library(BlazeIterative INTERFACE)

#==========================================
# Set BlazeIterative source files
#==========================================
file(GLOB SOURCE_FILES *${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)

# target_sources not necessary for header-only library, but makes it more
# IDE-friendly. Need to have it trigger only in "build" mode (not install)
# in order to avoid an error triggering (cmake v 3.9)
target_sources(BlazeIterative INTERFACE $<BUILD_INTERFACE:${SOURCE_FILES}>)
target_include_directories(BlazeIterative INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:include/>
)

#==========================================
# Import Blaze and "link" to BlazeIterative
#==========================================
find_package(blaze REQUIRED)
target_link_libraries(BlazeIterative INTERFACE blaze::blaze)
find_package(LAPACK REQUIRED)
target_link_libraries(BlazeIterative INTERFACE  ${LAPACK_LIBRARIES})

#==========================================
# Install library
#==========================================
install(TARGETS BlazeIterative EXPORT BlazeIterativeConfig)
install(
  DIRECTORY include/
  DESTINATION include/
  FILES_MATCHING PATTERN "*.hpp")
install(EXPORT BlazeIterativeConfig DESTINATION share/BlazeIterative/cmake
  NAMESPACE BlazeIterative::)

if(BUILD_TESTS)
  ENABLE_TESTING()
  add_subdirectory(tests)
endif()
