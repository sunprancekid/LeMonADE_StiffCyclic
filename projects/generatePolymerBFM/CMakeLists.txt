# add executables and compile
# generate polymer configuration and BFM file for simulation
add_executable(generatePolymerBFM generatePolymerBFM.cpp)
target_link_libraries(generatePolymerBFM LeMonADE)
