cmake_minimum_required(VERSION 3.0.0)
project(eigentools VERSION 0.1.0)

include(CTest)
enable_testing()

include_directories("C:/Users/735011/workspace/3rdparty/boost")
include_directories("C:/Users/735011/workspace/3rdparty/eigen3")

add_library(fastlm fastLm.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
