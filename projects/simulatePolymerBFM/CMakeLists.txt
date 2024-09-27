# add executables and compile
# real polymer simulation
add_executable(simulatePolymerBFM simulatePolymerBFM.cpp)
target_link_libraries(simulatePolymerBFM LeMonADE)

# ideal polymer simulation
add_executable(simulateIdealPolymerBFM simulateIdealPolymerBFM.cpp)
target_link_libraries(simulateIdealPolymerBFM LeMonADE)
