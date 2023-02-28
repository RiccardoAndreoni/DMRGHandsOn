cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(quantum_mc LANGUAGES CXX)

# set (CMAKE_CXX_FLAGS_DEBUG "-Wfatal-errors -Wparentheses -Wall -Wextra -Wpedantic -Wconversion -Werror -pedantic -std=c++17 -O0 -p -g -DDEBUG")
# message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")

# set (CMAKE_CXX_FLAGS "-Wfatal-errors -Wparentheses -Wall -Wextra -Wpedantic -Wconversion -Werror -pedantic -std=c++17 -O3")
set (CMAKE_CXX_FLAGS "-std=c++17 -O3")
message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")

include_directories(include/) # My libraries

add_library(
	common
	STATIC
    modules/block.cpp
    modules/sys.cpp
)

#================================================================================================================================
# Main programs

set(EXECUTABLE_OUTPUT_PATH ..)

add_executable(run_DMRG main/run_DMRG.cpp)
target_link_libraries(run_DMRG common) # pthread)

#================================================================================================================================
# Check programs
