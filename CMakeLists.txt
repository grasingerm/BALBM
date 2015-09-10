cmake_minimum_required(VERSION 2.8)

project(BALBM)

# set C++ standard
set(STDCXX c++14)
#set(STDCXX c++11)

# set flags for warnings, errors, and standards
set(CMAKE_BUILD_TYPE debug)
set(CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -g -std=c90 -pedantic -pedantic-errors")
set(CMAKE_C_FLAGS_RELEASE "-Wall -Wextra -std=c90 -pedantic -pedantic-errors -O3")
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g -std=${STDCXX} -pedantic -pedantic-errors")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -std=${STDCXX} -pedantic -pedantic-errors -O3")

# include and link directories
include_directories(include)

# dependencies
add_subdirectory(test)
