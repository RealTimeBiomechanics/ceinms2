add_executable(sandboxMuscleModel mainSandboxMuscleModel.cpp)
target_link_libraries(sandboxMuscleModel modeling project_options  ${Boost_LIBRARIES} )
target_include_directories(sandboxMuscleModel PRIVATE ${Boost_INCLUDE_DIRS})

add_executable(sandboxMuscleMapping mainSandboxMuscleMapping.cpp)
target_link_libraries(sandboxMuscleMapping PRIVATE modeling project_options  ${Boost_LIBRARIES} Eigen3::Eigen)
target_include_directories(sandboxMuscleMapping PRIVATE ${Boost_INCLUDE_DIRS})

add_executable(sandbox mainSandbox.cpp)
target_link_libraries(sandbox PRIVATE modeling project_options  ${Boost_LIBRARIES} Eigen3::Eigen)
target_include_directories(sandbox PRIVATE ${Boost_INCLUDE_DIRS})