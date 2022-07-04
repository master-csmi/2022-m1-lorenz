find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

# Find all the public headers
file(GLOB_RECURSE LORENZ_PUBLIC_HEADERS ${PROJECT_SOURCE_DIR}/src/cpp)

set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/src/cpp)
set(DOXYGEN_MAINPAGE_DIR ${PROJECT_SOURCE_DIR}/docs/doxygen)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/xml/index.xml)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/doxygenDoxyfile)

# Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

# Doxygen won't create this for us
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

# Only regenerate Doxygen when the Doxyfile or public headers change
add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
    COMMENT "Building doxygen..."
    DEPENDS ${CAT_CUTIFIER_PUBLIC_HEADERS}
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
    MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
    VERBATIM)

# Nice named target so we can run the job easily
add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)
configure_file(${SPHINX_SOURCE}/conf.py.in ${SPHINX_SOURCE}/conf.py @ONLY)
# Only regenerate Sphinx when:
# - Doxygen has rerun
# - Our doc files have been updated
# - The Sphinx config has been updated

add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
    COMMENT "Building sphinx..."
    COMMAND
    ${SPHINX_EXECUTABLE} -b html

    # # Tell Breathe where to find the Doxygen output
    # -Dbreathe_projects.lorenz=${DOXYGEN_OUTPUT_DIR}/xml
    ${SPHINX_SOURCE} ${SPHINX_BUILD}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS

    # Other docs files you want to track should go here (or in some variable)
    ${SPHINX_SOURCE}/index.rst
    # ${DOXYGEN_INDEX_FILE}
    MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py)

# Nice named target so we can run the job easily
add_custom_target(Sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})

# Add an install target to install the docs
include(GNUInstallDirs)
install(DIRECTORY ${SPHINX_BUILD}
    DESTINATION ${CMAKE_INSTALL_DOCDIR})

add_subdirectory(antora)

add_custom_target(docs ALL DEPENDS Doxygen Sphinx Antora)