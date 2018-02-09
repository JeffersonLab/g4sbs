#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(G4SBS)


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
# Update header file
add_custom_target( _gitinfo ALL COMMAND ${PROJECT_SOURCE_DIR}/pullgitinfo.py ${PROJECT_SOURCE_DIR} COMMENT "Creating gitinfo.hh" )


#ROOT path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
find_package(ROOT REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/src/dss2007 ${CMAKE_CURRENT_BINARY_DIR}/include ${PROJECT_SOURCE_DIR}/RTPC_2v4/include/ )
 
ROOT_GENERATE_DICTIONARY(
    G4SBSDict                         # path to dictionary to generate
    "include/G4SBSRunData.hh"         # list of classes to process
    "include/G4SBSTextFile.hh"
    LINKDEF "include/G4SBSLinkDef.h"  # ROOT linkDef file
)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/G4SBS*.cc ${PROJECT_SOURCE_DIR}/src/dss2007/*.cc ${PROJECT_SOURCE_DIR}/include/Pythia6_tree.C ${PROJECT_SOURCE_DIR}/RTPC_2v4/src/RTPC.cc ${PROJECT_SOURCE_DIR}/RTPC_2v4/src/Array*.cc ${PROJECT_SOURCE_DIR}/RTPC_2v4/src/RootIO.cc ${PROJECT_SOURCE_DIR}/RTPC_2v4/src/*Field*.cc ${PROJECT_SOURCE_DIR}/RTPC_2v4/src/BeamLine.cc )
file(GLOB headers ${CMAKE_CURRENT_BINARY_DIR}/include/gitinfo.hh ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/include/wiser_pion.h ${PROJECT_SOURCE_DIR}/src/cteq/*.h ${PROJECT_SOURCE_DIR}/src/dss2007/*.hh ${PROJECT_SOURCE_DIR}/include/Pythia6_tree.h ${PROJECT_SOURCE_DIR}/RTPC_2v4/include/RTPC.hh ${PROJECT_SOURCE_DIR}/RTPC_2v4/include/Array*.cc ${PROJECT_SOURCE_DIR}/RTPC_2v4/include/RootIO.hh ${PROJECT_SOURCE_DIR}/RTPC_2v4/include/*Field*.hh)

file(GLOB cteqsources ${PROJECT_SOURCE_DIR}/src/cteq/*.c  )
file(GLOB cteqheaders ${PROJECT_SOURCE_DIR}/src/cteq/*.h )

file(GLOB libsources ${PROJECT_SOURCE_DIR}/src/G4SBSRunData.cc ${PROJECT_SOURCE_DIR}/src/G4SBSTextFile.cc)
file(GLOB libheaders ${CMAKE_CURRENT_BINARY_DIR}/include/gitinfo.hh ${PROJECT_SOURCE_DIR}/include/RTPC.hh )

add_subdirectory(RTPC_2v4)

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

include(${Geant4_USE_FILE})
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(sbscteq STATIC ${cteqsources} ${cteqheaders})
SET_TARGET_PROPERTIES(sbscteq PROPERTIES COMPILE_FLAGS "-DCTEQ_TBL_PATH=\\\"${PROJECT_SOURCE_DIR}/cteq-tbls\\\"")

add_executable(g4sbs g4sbs.cc ${sources} ${headers} G4SBSDict.cxx)
add_library(g4sbsroot SHARED ${libsources} ${libheaders} G4SBSDict.cxx)

target_link_libraries(g4sbs ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} sbscteq )
target_link_libraries(g4sbsroot ${ROOT_LIBRARIES} )

option(WITH_G4SBS_GDML "Build g4sbs with GDML output" OFF)
if( WITH_G4SBS_GDML)
    set_target_properties(g4sbs 
    PROPERTIES COMPILE_FLAGS "-DG4SBS_USE_GDML")
endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build g4sbs. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(G4SBS_SCRIPTS
  vis.mac
  run_example.mac
  scripts/gen_1.5GeV2.mac
  scripts/gen_4GeV2.mac
  scripts/gen_7GeV2.mac
  scripts/gen_10GeV2.mac
  scripts/gmn_3.5GeV2.mac
  scripts/gmn_4.5GeV2.mac
  scripts/gmn_5.7GeV2.mac
  scripts/gmn_8.1GeV2.mac
  scripts/gmn_10.2GeV2.mac
  scripts/gmn_12GeV2.mac
  scripts/gmn_13.5GeV2.mac
  scripts/gep_5GeV2.mac
  scripts/gep_8GeV2.mac
  scripts/gep_12GeV2.mac
  #scripts/sidis_template.mac
  #scripts/make_SIDIS_scripts.C
  scripts/vis_sidis.mac
  scripts/sidis_background.mac
  scripts/vis_gep.mac
  scripts/vis_gmn.mac
  scripts/preinit_ckov_scint.mac
  scripts/preinit_ckov_noscint.mac
  scripts/preinit_nockov_scint.mac
  scripts/preinit_nockov_noscint.mac
  scripts/sidis_example_script.mac
  root_macros/G4SBS_optics_fit.C 
  root_macros/RICH_analysis.C    
  root_macros/g4sbs_tree.C       
  root_macros/g4sbs_tree.h       
  root_macros/rate_table.C
  root_macros/SIDIS_reduced_data_file.C
  root_macros/SIDIS_DST.C
  root_macros/SIDIS_DST.h
  root_macros/rate_table_DST.C
  root_macros/Pythia6_minbias_gen.C
  database/ecal_map.txt
  scripts/preinit_ckov_scint.mac
  scripts/preinit_ckov_noscint.mac
  scripts/preinit_nockov_scint.mac
  scripts/preinit_nockov_noscint.mac
  )

foreach(_script ${G4SBS_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

configure_file( ${PROJECT_SOURCE_DIR}/rootlogon.C ${PROJECT_BINARY_DIR}/rootlogon.C )

set(DSS2007_GRIDS HLO.GRID HNLO.GRID KALO.GRID KNLO.GRID PILO.GRID PROLO.GRID PLO.GRID PNLO.GRID)

foreach(_grid ${DSS2007_GRIDS})
  configure_file( ${PROJECT_SOURCE_DIR}/src/dss2007/${_grid} ${PROJECT_BINARY_DIR}/${_grid} COPYONLY )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(G4SBS DEPENDS g4sbs)

add_dependencies(g4sbs _gitinfo)
add_dependencies(g4sbsroot _gitinfo)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS g4sbs DESTINATION bin)


