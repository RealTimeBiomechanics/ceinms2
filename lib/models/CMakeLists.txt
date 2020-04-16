set(models_HEADERS  ceinms2/Lloyd2019Muscle.h
                    ceinms2/Lloyd2019Activation.h
					ceinms2/WDBsolver.h
					)
set(models_TEMPLATE_IMPLEMENTATION
					)
set(models_SOURCES
                    models.cpp)
set_source_files_properties( ${models_TEMPLATE_IMPLEMENTATION} PROPERTIES HEADER_FILE_ONLY TRUE)
source_group("Header Files" FILES ${models_HEADERS})
source_group("Source Files" FILES ${models_TEMPLATE_IMPLEMENTATION})


# Save path to include directory for "internal" usage
set(models_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "Include directory for models library")

add_library(models STATIC ${models_HEADERS} ${models_TEMPLATE_IMPLEMENTATION} ${models_SOURCES})
set_target_properties(models PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(models PUBLIC common)
target_include_directories(models PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS models
	EXPORT CEINMSTargets
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install
install(FILES ${models_HEADERS} ${models_TEMPLATE_IMPLEMENTATION} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ceinms2 COMPONENT Development)
