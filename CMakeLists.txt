cmake_minimum_required(VERSION 3.14)
project(ufo VERSION 0 LANGUAGES CXX)
enable_testing()
include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message("Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

find_package(TBB REQUIRED)
find_package(Matplot++ REQUIRED)
find_package(biu REQUIRED)
find_package(Threads REQUIRED)
message("biu dir: ${biu_DIR}")

add_executable(ufo src/fold.cpp src/unfold.cpp src/project-to-mode.cpp
  src/plot.cpp src/raman.cpp src/project-to-atom.cpp src/main.cpp)
target_include_directories(ufo PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(ufo PRIVATE TBB::tbb Matplot++::matplot biu::biu)
target_compile_features(ufo PRIVATE cxx_std_23)
target_compile_definitions(ufo PRIVATE BIU_LOGGER_SOURCE_ROOT="${CMAKE_SOURCE_DIR}")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_definitions(ufo PRIVATE BIU_LOGGER_DEBUG)
  target_compile_options(ufo PRIVATE -ftemplate-backtrace-limit=0)
endif()

install(TARGETS ufo RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

get_property(ImportedTargets DIRECTORY "${CMAKE_SOURCE_DIR}" PROPERTY IMPORTED_TARGETS)
message("Imported targets: ${ImportedTargets}")
message("List of compile features: ${CMAKE_CXX_COMPILE_FEATURES}")
message("CMake build type: ${CMAKE_BUILD_TYPE}")

include(CTest)

add_test(NAME fold WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/fold COMMAND ufo fold config.yaml)
add_test(NAME unfold WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/unfold COMMAND ufo unfold config.yaml)
add_test(NAME zpp-to-yaml WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/zpp-to-yaml
  COMMAND sh -c "${CMAKE_BINARY_DIR}/ufo zpp-to-yaml < data.zpp > data.yaml")
add_test(NAME project-to-mode WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/project-to-mode
  COMMAND ufo project-to-mode config.yaml)
add_test(NAME project-to-atom WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/project-to-atom
  COMMAND ufo project-to-atom config.yaml)
add_test(NAME raman-create-displacement WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/raman-create-displacement
  COMMAND ufo raman-create-displacement config.yaml)
add_test(NAME raman-extract WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/raman-extract
  COMMAND ufo raman-extract .)
add_test(NAME raman-apply-contribution WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/raman-apply-contribution
  COMMAND ufo raman-apply-contribution config.yaml)
