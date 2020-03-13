cmake_minimum_required(VERSION 3.12)

project(mibrr
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(eigen3)
find_package(Boost)

include_directories(${Boost_INCLUDE_DIRS})

add_executable(mibrr src/main.cpp src/MibrrData.h src/MibrrGibbs.h src/MibrrSamplers.h
					 src/MibrrData.cpp src/MibrrGibbs.cpp
					 src/MibrrSamplers.cpp)

target_link_libraries(mibrr Eigen3::Eigen ${BOOST_LIBRARIES})