cmake_minimum_required(VERSION 3.0.2)
project(toppra_optimizer_ros)

add_compile_options(-std=c++11)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  geometry_msgs
  std_msgs
  trajectory_msgs
  genmsg
  actionlib_msgs
  actionlib
)



find_package(Boost REQUIRED COMPONENTS system)
find_package(Eigen3)

find_package(toppra)
find_package(NumCpp)


add_action_files(
  DIRECTORY action
  FILES Optimize.action
)

generate_messages(
  DEPENDENCIES actionlib_msgs trajectory_msgs std_msgs geometry_msgs
)

catkin_package(
  CATKIN_DEPENDS actionlib_msgs
)

catkin_package(
  LIBRARIES toppra
)


option(PYTHON_BINDINGS "Build bindings for Python" ON)
set(PYTHON_BINDINGS OFF)

set(PYTHON_VERSION 3.8 CACHE STRING "Build bindings for Python version")

option(BUILD_TESTS "Build Gtests" ON)
set(BUILD_TESTS OFF)

# add_subdirectory(toppra/cpp)

include_directories(${EIGEN_INCLUDE_DIRS})

# add_action_files(DIRECTORY action FILES ACTION_FOR_OPTIMIZATION.action) TODO


include_directories(

  ${catkin_INCLUDE_DIRS}
)


add_executable(toppra_optimizer_server src/toppra_optimizer_server.cpp)

target_link_libraries(
  toppra_optimizer_server
  PUBLIC toppra::toppra
  ${catkin_LIBRARIES}
)
add_dependencies(
  toppra_optimizer_server
  ${toppra_optimizer_ros_EXPORTED_TARGETS}
)