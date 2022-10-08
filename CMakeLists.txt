cmake_minimum_required (VERSION 3.8)

project ("Heat Equationn")

add_executable (Solver "Heat Equationn.cpp" "Heat Equationn.h" "Function.cpp" "Function.h" "Parser.h" "Parser.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Solver PROPERTY CXX_STANDARD 20)
endif()
