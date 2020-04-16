# Author: Elena Ceseracciu

# List source files
set(models_HEADERS  ceinms2/Lloyd2019Muscle.h
                    ceinms2/Lloyd2019Activation.h)
set(models_TEMPLATE_IMPLEMENTATION
                    )

set_source_files_properties( ${models_TEMPLATE_IMPLEMENTATION} PROPERTIES HEADER_FILE_ONLY TRUE)

# Save path to include directory for "internal" usage
set(models_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "Include directory for models library")

if (NOT CMAKE_VERSION VERSION_LESS 3.0.0)
     add_library(models INTERFACE)
     install(TARGETS models
        EXPORT CEINMSTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# Install
install(FILES ${models_HEADERS} ${models_TEMPLATE_IMPLEMENTATION} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ceinms2 COMPONENT Development)
