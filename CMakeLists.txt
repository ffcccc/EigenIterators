cmake_minimum_required(VERSION 3.0.0)
project(eigentools VERSION 0.1.0)

link_directories( build )

set ( PROJECT_INC_BOOST C:/Users/735011/workspace/3rdparty/boost )
set ( PROJECT_INC_EIGEN C:/Users/735011/workspace/3rdparty/eigen3 )
set ( PROJECT_INC_UTILS C:/Users/735011/workspace/KCL/EigenUtils )
set ( PROJECT_INC_CSV C:/Users/735011/workspace/3rdparty/fast-csv-parser-master )
set ( PROJECT_INCS C:/Users/735011/workspace/3rdparty )

include_directories(${PROJECT_INCS})
include_directories(${PROJECT_INC_BOOST} ${PROJECT_INC_EIGEN} ${PROJECT_INC_UTILS} ${PROJECT_INC_CSV} ${PROJECT_INCS} )

add_executable(eigentools main.cpp)