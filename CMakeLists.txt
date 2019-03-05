cmake_minimum_required(VERSION 3.13)
project(ThreesACrowd)

set(CMAKE_CXX_STANDARD 14)

# https://github.com/jgaeddert/liquid-dsp/issues/66
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -msse3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS-RELEASE} -msse3")

include_directories(src, glm)

add_executable(ThreesACrowd
        .idea/encodings.xml
        .idea/misc.xml
        .idea/modules.xml
        .idea/ThreesACrowd.iml
        .idea/vcs.xml
        .idea/workspace.xml
        src/Grid.cpp
        src/Grid.h
        src/main.cpp
        README.md src/Person.cpp src/Person.h)
