add_executable(sandboxMuscleModel mainSandboxMuscleModel.cpp)

target_link_libraries(sandboxMuscleModel modeling project_options  ${Boost_LIBRARIES} )
target_include_directories(sandboxMuscleModel PRIVATE ${Boost_INCLUDE_DIRS})