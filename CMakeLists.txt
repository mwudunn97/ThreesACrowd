cmake_minimum_required(VERSION 3.13)
project(ThreesACrowd)

set(CMAKE_CXX_STANDARD 14)

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
        README.md)
