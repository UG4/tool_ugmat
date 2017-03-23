################################################
# created by Sebastian Reiter
# s.b.reiter@gmail.com
################################################

cmake_minimum_required(VERSION 2.6)
project(matutil)

set(sources	matutil.cpp)

add_executable(matutil ${sources})
install(TARGETS matutil RUNTIME DESTINATION "bin")
