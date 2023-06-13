cmake_minimum_required(VERSION 3.12)

project(AssemblyTesting)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")

add_executable(tester tester/main.cpp tester/assembly.cpp tester/assembly.h
        tester/sequencing.cpp tester/sequencing.h tester/utility.cpp tester/utility.h)

add_executable(solver solver/main.cpp solver/simulated_annealing.cpp
        solver/simulated_annealing.h solver/utility.cpp solver/utility.h)

add_executable(meter meter/main.cpp meter/levenshtein.cpp meter/levenshtein.h meter/ssw.c meter/ssw.h
        meter/ssw_cpp.cpp meter/ssw_cpp.h)
