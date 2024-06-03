cmake_minimum_required (VERSION 2.8)

if (NOT DEFINED LEMONADE_DIR)
    message("LEMONADE_DIR is not provided.LEMONADE_DIR is the installation path of the main repo. If build fails, use -DLEMONADE_DIR=/path/to/LeMonADE/")
endif()

#specify the paths to find the include files and the library (.cpp files)
set (LEMONADE_INCLUDE_DIR ${LEMONADE_DIR}include/)
set (LEMONADE_LIBRARY_DIR ${LEMONADE_DIR}lib/)

message( "LeMonADE_installPath : " ${LEMONADE_DIR} )
message( "LeMonADE_includePath : " ${LEMONADE_INCLUDE_DIR} )
message( "LeMonADE_libraryPath : " ${LEMONADE_LIBRARY_DIR} )

include_directories (${LEMONADE_INCLUDE_DIR})
link_directories (${LEMONADE_LIBRARY_DIR})
set (EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

#define flags for the gcc compiler 
set (CMAKE_CXX_FLAGS " -O3 -msse2 -mssse3 -std=c++11 -fexpensive-optimizations -std=c++11 -Wno-error=narrowing")
set (CMAKE_C_FLAGS "-O3 -msse2 -mssse3 -std=c++11-fexpensive-optimizations -std=c++11 -Wno-error=narrowing")

#add the relevant subdirectories for doing fancy stuff 
add_subdirectory(projects)
# add_subdirectory(analyzer)
# add_subdirectory(updater)
# add_subdirectory(feature)


if (NOT DEFINED LEMONADE_TESTS )
    message("If you want to compile the tests use -DLEMONADE_TESTS=ON")
endif()
# option(LEMONADE_TESTS "Build the test" OFF)
if(LEMONADE_TESTS)
    message("Add the ./test directory and compile the tests")
    add_subdirectory(test)
endif(LEMONADE_TESTS)
