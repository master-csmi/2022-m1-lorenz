set(APPLICATION_DOC_DIR "${CMAKE_INSTALL_DOCDIR}")

set(ANTORA_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/public")
message("${ANTORA_OUTPUT_DIR}")
find_program (NPM_EXECUTABLE NAMES npm)

# Doxygen won't create this for us
file(MAKE_DIRECTORY ${ANTORA_OUTPUT_DIR})

if ( NPM_EXECUTABLE )
    add_custom_command(OUTPUT "${ANTORA_OUTPUT_DIR}" ALL
        COMMENT "Building documentation... "
        COMMAND "${npm}" install 
        COMMAND "${npm}" run antora 
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )

    if ( NOT TARGET Antora )
        add_custom_target(Antora ALL DEPENDS "${ANTORA_OUTPUT_DIR}")
    else()
        add_dependencies(Antora "${ANTORA_OUTPUT_DIR}")
    endif()

    install(DIRECTORY "${ANTORA_OUTPUT_DIR}/" DESTINATION ${APPLICATION_DOC_DIR} OPTIONAL)
else()
    message(STATUS "npm is not installed, can use antora ")
endif()