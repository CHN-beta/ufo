cmake_minimum_required(VERSION 3.14)
project(ufo VERSION 0.0.0 LANGUAGES CXX)
enable_testing()
include(GNUInstallDirs)

set_property(GLOBAL PROPERTY CXX_STANDARD 23)
set_property(GLOBAL PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(GLOBAL PROPERTY CXX_EXTENSIONS OFF)

find_package(yaml-cpp REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(fmt REQUIRED)
find_package(concurrencpp REQUIRED)
find_package(HighFive REQUIRED)
find_package(TBB REQUIRED)
find_package(Matplot++ REQUIRED)
find_path(ZPP_BITS_INCLUDE_DIR zpp_bits.h REQUIRED)

add_executable(ufo src/ufo.cpp)
target_include_directories(ufo PRIVATE ${PROJECT_SOURCE_DIR}/include ${ZPP_BITS_INCLUDE_DIR})
target_link_libraries(ufo PRIVATE
  yaml-cpp Eigen3::Eigen fmt::fmt concurrencpp::concurrencpp HighFive_HighFive TBB::tbb Matplot++::matplot)

get_property(ImportedTargets DIRECTORY "${CMAKE_SOURCE_DIR}" PROPERTY IMPORTED_TARGETS)
message("Imported targets: ${ImportedTargets}")
message("List of compile features: ${CMAKE_CXX_COMPILE_FEATURES}")
