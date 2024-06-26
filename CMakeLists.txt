cmake_minimum_required(VERSION 3.14)
project(ufo VERSION 0.0.0 LANGUAGES CXX)
enable_testing()
include(GNUInstallDirs)

set_property(GLOBAL PROPERTY CXX_STANDARD 23)
set_property(GLOBAL PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(GLOBAL PROPERTY CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message("Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

find_package(yaml-cpp REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(fmt REQUIRED)
find_package(concurrencpp REQUIRED)
find_package(HighFive REQUIRED)
find_package(TBB REQUIRED)
find_package(glad REQUIRED)
# find_package(Matplot++ REQUIRED COMPONENTS matplot_opengl)
find_package(Matplot++ REQUIRED)
find_path(ZPP_BITS_INCLUDE_DIR zpp_bits.h REQUIRED)

add_executable(ufo src/solver.cpp src/fold.cpp src/unfold.cpp src/plot.cpp src/main.cpp)
target_include_directories(ufo PRIVATE ${PROJECT_SOURCE_DIR}/include ${ZPP_BITS_INCLUDE_DIR})
target_link_libraries(ufo PRIVATE
  yaml-cpp Eigen3::Eigen fmt::fmt concurrencpp::concurrencpp HighFive_HighFive TBB::tbb Matplot++::matplot
  Matplot++::matplot_opengl)
# target_compile_definitions(ufo PRIVATE $<$<CONFIG:Debug>:_GLIBCXX_DEBUG>)

get_property(ImportedTargets DIRECTORY "${CMAKE_SOURCE_DIR}" PROPERTY IMPORTED_TARGETS)
message("Imported targets: ${ImportedTargets}")
message("List of compile features: ${CMAKE_CXX_COMPILE_FEATURES}")
