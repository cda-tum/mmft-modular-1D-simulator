# enable_testing()
# include(GoogleTest)

# add_executable(dropletTest)
# target_link_libraries(${TARGET_NAME} PUBLIC gtest gtest_main)
# target_link_libraries(${TARGET_NAME} PRIVATE droplet)

add_subdirectory(architecture)
add_subdirectory(nodalAnalysis)
add_subdirectory(simulation)
add_subdirectory(gradientgenerator)
add_subdirectory(organ)
add_subdirectory(droplet)

# gtest_discover_tests(${TARGET_NAME})