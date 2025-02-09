# Set the minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(VPTComm VERSION 1.0)

# Set the path for the include directory
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

# Set the MPI compiler to mpicc
set(CMAKE_C_COMPILER mpicc.mpich)

# Add the include directory
include_directories(${INCLUDE_DIR})

# Collect all .c files in the src directory
file(GLOB SOURCES "src/*.c")

# Create the static library libvptcomm.a
add_library(vptcomm STATIC ${SOURCES})

# Link the math library to the static library
target_link_libraries(vptcomm PRIVATE m)

# Specify the installation directory for the library
install(TARGETS vptcomm DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)

# Add all test source files (assuming they are in the 'tests' directory)
file(GLOB TEST_SOURCES "tests/*.c")

# Define each test executable and link it to the library
foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)  # Extract base name of the test file
    add_executable(${test_name} ${test_src})
    target_link_libraries(${test_name} vptcomm m)  # Link each test with vptcomm library and math library
    target_include_directories(${test_name} PRIVATE ${INCLUDE_DIR})  # Add include directory for tests
    # Register the test with CTest using mpirun -np 8
    add_test(NAME ${test_name} 
        COMMAND mpiexec.mpich -np 8 $<TARGET_FILE:${test_name}>)
endforeach()

# Enable testing with CTest
enable_testing()

# Add a custom target to run all tests after build
add_custom_target(run_tests
    COMMAND ${CMAKE_CTEST_COMMAND}
    DEPENDS ${TEST_SOURCES}
)

# Optionally, add a 'make install' target for the package
install(TARGETS vptcomm DESTINATION lib)

