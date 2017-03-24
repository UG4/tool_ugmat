# This file is part of ugmat, a program for analysing and comparing matrices
#
# Copyright (C) 2017 Sebastian Reiter, G-CSC Frankfurt <sreiter@gcsc.uni-frankfurt.de>
#
# This Source Code Form is subject to the terms of the Mozilla
# Public License v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 2.6)
project(ugmat)

set(sources	ugmat.cpp)

add_executable(ugmat ${sources})
install(TARGETS ugmat RUNTIME DESTINATION "bin")
