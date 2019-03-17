cmake_minimum_required(VERSION 3.12)
project(ThreesACrowd)

set(CMAKE_CXX_STANDARD 14)

# https://github.com/jgaeddert/liquid-dsp/issues/66
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -msse3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS-RELEASE} -msse3 -O3")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

# Add source files
file(GLOB_RECURSE SOURCE_FILES
        ${CMAKE_SOURCE_DIR}/src/*.c
        ${CMAKE_SOURCE_DIR}/src/*.cpp)

# Add header files
file(GLOB_RECURSE HEADER_FILES
        ${CMAKE_SOURCE_DIR}/src/*.h
        ${CMAKE_SOURCE_DIR}/src/*.hpp)

include_directories(src glm nlohmann)

add_executable(ThreesACrowd
        .idea/encodings.xml
        .idea/misc.xml
        .idea/modules.xml
        .idea/ThreesACrowd.iml
        .idea/vcs.xml
        .idea/workspace.xml
        README.md
        ${HEADER_FILES}
        ${SOURCE_FILES} src/pointDisplay.h)

target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

