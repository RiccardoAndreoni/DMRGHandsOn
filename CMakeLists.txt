cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(DMRG LANGUAGES CXX)

# set (CMAKE_CXX_FLAGS_DEBUG "-Wfatal-errors -Wparentheses -Wall -Wextra -Wpedantic -Wconversion -Werror -pedantic -std=c++17 -O0 -p -g -DDEBUG")
# message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")

set (CMAKE_CXX_FLAGS "-Wfatal-errors -Wparentheses -Wall -Wextra -Wpedantic -Wconversion -Werror -pedantic -std=c++17 -O3")
# set (CMAKE_CXX_FLAGS "-std=c++17 -O3")
message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")

include_directories(include/) # My libraries

add_library(
	common
	STATIC
    modules/block.cpp
    modules/sys.cpp
    modules/model.cpp
    modules/DMRG.cpp
    modules/service.cpp
    )

#================================================================================================================================
# Main programs

set(EXECUTABLE_OUTPUT_PATH)

# find_package(GSL REQUIRED)
# add_executable(run_DMRG main/run_DMRG.cpp)
# target_link_libraries(run_DMRG common GSL::gsl GSL::gslcblas) # pthread)

#================================================================================================================================
# Check programs

# find_package(GSL REQUIRED)
# add_executable(test_tp test/test_tp.cpp)
# target_link_libraries(test_tp common GSL::gsl GSL::gslcblas) # pthread)

# find_package(GSL REQUIRED)
# add_executable(test_reshape test/test_reshape.cpp)
# target_link_libraries(test_reshape common GSL::gsl GSL::gslcblas) # pthread)

# find_package(GSL REQUIRED)
# add_executable(test_lanczos1 test/test_lanczos1.cpp)
# target_link_libraries(test_lanczos1 common GSL::gsl GSL::gslcblas) # pthread)

# find_package(GSL REQUIRED)
# add_executable(test_model test/test_model.cpp)
# target_link_libraries(test_model common GSL::gsl GSL::gslcblas) # pthread)

# find_package(GSL REQUIRED)
# add_executable(test_add_site test/test_add_site.cpp)
# target_link_libraries(test_add_site common GSL::gsl GSL::gslcblas) # pthread)

# find_package(GSL REQUIRED)
# add_executable(test_convert test/test_convert.cpp)
# target_link_libraries(test_convert common GSL::gsl GSL::gslcblas) # pthread)

# find_package(GSL REQUIRED)
# add_executable(test_SVD test/test_SVD.cpp)
# target_link_libraries(test_SVD common GSL::gsl GSL::gslcblas) # pSVD

# find_package(GSL REQUIRED)
# add_executable(test_truncation test/test_truncation.cpp)
# target_link_libraries(test_truncation common GSL::gsl GSL::gslcblas) # pthread)

# find_package(GSL REQUIRED)
# add_executable(test_pointers1 test/test_pointers1.cpp)
# target_link_libraries(test_pointers1 common GSL::gsl GSL::gslcblas) # pthread)

find_package(GSL REQUIRED)
add_executable(test_renormalize test/test_renormalize.cpp)
target_link_libraries(test_renormalize common GSL::gsl GSL::gslcblas) # pthread)

# find_package(GSL REQUIRED)
# add_executable(test_lanczos2 test/test_lanczos2.cpp)
# target_link_libraries(test_lanczos2 common GSL::gsl GSL::gslcblas) # pthread)


# find_package(GSL REQUIRED)
# add_executable(test_infinite test/test_infinite.cpp)
# target_link_libraries(test_infinite common GSL::gsl GSL::gslcblas) # pthread)