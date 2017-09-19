cmake_minimum_required(VERSION 3.0)
project(kseqcpp)

set(CMAKE_CXX_STANDARD 11)

find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
    add_definitions(-DHAVE_ZLIB)
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

find_package(BZip2 REQUIRED)
if (BZIP2_FOUND)
    add_definitions(-DHAVE_BZIP2)
    include_directories(${BZIP2_INCLUDE_DIRS})
endif()

set(SOURCE_FILES test.cpp kseq.hpp)

add_executable(kseq ${SOURCE_FILES})

target_link_libraries (kseq ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES})
