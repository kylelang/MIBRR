cmake_minimum_required(VERSION 3.12)

project(mibrr
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(eigen3 REQUIRED)
find_package(Boost)
find_package(nlohmann_json REQUIRED)
find_package(spdlog REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

set(NLOHMANN_JSON_LIBRARIES "nlohmann_json::nlohmann_json")
set(SPDLOG_LIBRARIES "spdlog::spdlog")

add_executable(mibrr src/main.cpp src/MibrrData.h src/MibrrGibbs.h src/MibrrSamplers.h
					 src/MibrrData.cpp src/MibrrGibbs.cpp src/MibrrRunGibbs.cpp
					 src/MibrrSamplers.cpp)

target_link_libraries(mibrr Eigen3::Eigen
							${BOOST_LIBRARIES}
							${NLOHMANN_JSON_LIBRARIES}
							${SPDLOG_LIBRARIES})