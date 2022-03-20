cmake_minimum_required(VERSION 3.13)
project(Learned-Compressed-Rank-Select-TALG22)

set(CMAKE_CXX_STANDARD 20)

configure_file(
        include/ds2i/ds2i_config.hpp.in
        include/ds2i/ds2i_config.hpp
        ESCAPE_QUOTES)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -pthread -fopenmp")

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.45.0 COMPONENTS program_options)

include_directories(include lib/la_vector/include include/ds2i lib)

add_subdirectory(lib/succinct)
add_subdirectory(lib/la_vector)

if (NOT Boost_FOUND)
    message(FATAL_ERROR "Fatal error: Boost required.")
endif (NOT Boost_FOUND)

include_directories(${Boost_INCLUDE_DIRS})

add_executable(my_benchmark my_benchmark.cpp)
target_link_libraries(my_benchmark sdsl divsufsort divsufsort64 ${Boost_LIBRARIES})

add_executable(time_to_build time_to_build.cpp)
target_link_libraries(time_to_build sdsl divsufsort divsufsort64 ${Boost_LIBRARIES})

add_subdirectory(tests)