# compile executable for generating a polymer configuration according to certain parameters
add_subdirectory(generatePolymerBFM)
# compule executable for simulation a generalized polymer configuration
add_subdirectory(simulatePolymerBFM)
# create two concatenated rings with excluded volume 
add_subdirectory(createTwoConcatenatedRings)
