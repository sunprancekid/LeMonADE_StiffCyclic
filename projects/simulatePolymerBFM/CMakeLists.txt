# add executables and compile
# real polymer simulation
add_executable(simulateRealPolymerBFM simulateRealPolymerBFM.cpp)
target_link_libraries(simulateRealPolymerBFM LeMonADE)

# ideal polymer simulation
add_executable(simulateIdealPolymerBFM simulateIdealPolymerBFM.cpp)
target_link_libraries(simulateIdealPolymerBFM LeMonADE)
