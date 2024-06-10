# add executables and compile
# ideal chain
# add_executable(linearChainIdeal lemonade_ISC_main.cpp)
# target_link_libraries(linearChainIdeal LeMonADE)

# real chain
add_executable(linearChainReal_forceExtension lemonade_RSC_main.cpp)
target_link_libraries(linearChainReal_forceExtension LeMonADE)
