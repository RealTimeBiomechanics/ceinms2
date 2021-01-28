set(modeling_HEADERS ceinms2/Lloyd2003Muscle.h
                     ceinms2/ExponentialActivation.h
					 ceinms2/WDBsolver.h
					 ceinms2/NMSmodel.h
					)
set(modeling_TEMPLATE_IMPLEMENTATION
					)
set(modeling_SOURCES
                    modeling.cpp)
set_source_files_properties( ${modeling_TEMPLATE_IMPLEMENTATION} PROPERTIES HEADER_FILE_ONLY TRUE)
source_group("Header Files" FILES ${modeling_HEADERS})
source_group("Source Files" FILES ${modeling_TEMPLATE_IMPLEMENTATION})


# Save path to include directory for "internal" usage
set(modeling_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "Include directory for modeling library")

add_library(modeling STATIC ${modeling_HEADERS} ${modeling_TEMPLATE_IMPLEMENTATION} ${modeling_SOURCES})
set_target_properties(modeling PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(modeling PUBLIC common)
target_include_directories(modeling PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> ${Boost_INCLUDE_DIRS})

install(TARGETS modeling
	EXPORT CEINMSTargets
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install
install(FILES ${modeling_HEADERS} ${modeling_TEMPLATE_IMPLEMENTATION} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ceinms2 COMPONENT Development)
