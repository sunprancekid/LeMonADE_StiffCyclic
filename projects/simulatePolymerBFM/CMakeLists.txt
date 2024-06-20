# add executables and compile
# generate polymer configuration and BFM file for simulation
add_executable(simulatePolymerBFM simulatePolymerBFM.cpp)
target_link_libraries(simulatePolymerBFM LeMonADE)
