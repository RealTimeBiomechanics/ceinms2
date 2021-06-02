set(common_HEADERS 	ceinms2/Curve.h
					ceinms2/CircularVector.h
					ceinms2/DataTable.h
					ceinms2/Types.h
					ceinms2/Mapper.h)
set(common_TEMPLATE_IMPLEMENTATION
					ceinms2/Curve.cpp
					ceinms2/CircularVector.cpp
					ceinms2/DataTable.cpp)
set(common_SOURCES
                    common.cpp Mapper.cpp)

set_source_files_properties(${common_TEMPLATE_IMPLEMENTATION} PROPERTIES HEADER_FILE_ONLY TRUE)
source_group("Header Files" FILES ${common_HEADERS})
source_group("Source Files" FILES ${common_TEMPLATE_IMPLEMENTATION})


# Save path to include directory for "internal" usage
set(common_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "Include directory for common library")

add_library(common STATIC ${common_HEADERS} ${common_TEMPLATE_IMPLEMENTATION} ${common_SOURCES})
#target_sources(common INTERFACE ${common_HEADERS} ${common_TEMPLATE_IMPLEMENTATION})
set_target_properties(common PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(common PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS common
	EXPORT CEINMSTargets
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install
install(FILES ${common_HEADERS} ${common_TEMPLATE_IMPLEMENTATION} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ceinms2 COMPONENT Development)
