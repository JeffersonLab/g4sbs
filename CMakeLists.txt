#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(G4SBS)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

#ROOT path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
find_package(ROOT REQUIRED)

find_package(CLHEP)

include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/src/ )

if(CLHEP_FOUND)
    include_directories(${CLHEP_INCLUDE_DIRS})
    link_directories(${CLHEP_LIBRARY_DIRS})
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/G4SBS*.cc )
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/src/cteq/*.h )

file(GLOB cteqsources ${PROJECT_SOURCE_DIR}/src/cteq/*.c  )
file(GLOB cteqheaders ${PROJECT_SOURCE_DIR}/src/cteq/*.h )

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(sbscteq STATIC ${cteqsources} ${cteqheaders})
SET_TARGET_PROPERTIES(sbscteq PROPERTIES COMPILE_FLAGS "-DCTEQ_TBL_PATH=\\\"${PROJECT_SOURCE_DIR}/cteq-tbls\\\"")

add_executable(g4sbs g4sbs.cc ${sources} ${headers})

target_link_libraries(g4sbs ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} sbscteq )
if(CLHEP_FOUND)
    #    target_link_libraries(g4sbs ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} sbscteq CLHEP )
    target_link_libraries(g4sbs ${CLHEP_LIBRARIES} )
else()
    #    target_link_libraries(g4sbs ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} sbscteq )
endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build g4sbs. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(G4SBS_SCRIPTS
  vis.mac
  )

foreach(_script ${G4SBS_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(G4SBS DEPENDS g4sbs)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS g4sbs DESTINATION bin)


