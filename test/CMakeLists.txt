set(testing_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "Include directory for testing")		
							
add_executable(testNMSmodel testNMSmodel.cpp)
target_include_directories(testNMSmodel PRIVATE ${testing_INCLUDE_DIRS})
target_link_libraries(testNMSmodel modeling project_options project_warnings )
add_test(TestNMSmodel testNMSmodel)

add_executable(testMuscleModels testMuscleModels.cpp)
target_include_directories(testMuscleModels PRIVATE ${testing_INCLUDE_DIRS})
target_link_libraries(testMuscleModels modeling project_options project_warnings)
add_test(TestMuscleModels testMuscleModels)

add_executable(testDof testDof.cpp)
target_include_directories(testDof PRIVATE ${testing_INCLUDE_DIRS})
target_link_libraries(testDof modeling project_options project_warnings)
add_test(TestDof testDof)

add_executable(testLinearActuator testLinearActuator.cpp)
target_include_directories(testLinearActuator PRIVATE ${testing_INCLUDE_DIRS})
target_link_libraries(testLinearActuator modeling project_options project_warnings)
add_test(TestLinearActuator testLinearActuator)

add_executable(testElectromechanicalDelay testElectromechanicalDelay.cpp)
target_include_directories(testElectromechanicalDelay PRIVATE ${testing_INCLUDE_DIRS})
target_link_libraries(testElectromechanicalDelay modeling project_options project_warnings)
add_test(TestElectromechanicalDelay testElectromechanicalDelay)

add_executable(testBoostAutodiff testBoostAutodiff.cpp)
target_link_libraries(testBoostAutodiff modeling project_options project_warnings)
target_include_directories(testBoostAutodiff PRIVATE ${Boost_INCLUDE_DIRS} ${testing_INCLUDE_DIRS})
add_test(TestBoostAutodiff testBoostAutodiff)
