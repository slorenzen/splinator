project(splinator)

add_executable(splinator timeseries.cpp main.cpp)

target_link_libraries(splinator gsl gslcblas)