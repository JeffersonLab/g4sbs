#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(G4SBS)

# Setup GNUInstallDirs for installation:
#----------------------------------------------------------------------------
include(GNUInstallDirs)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# The following lines are required so that the executable can dynamically locate the library dependencies independent of 
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
ENDIF("${isSystemDir}" STREQUAL "-1")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------

#
# Force use of c++11 standard
#
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#

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
#include_directories(${PROJECT_SOURCE_DIR}/include)

#check for absence of GEANT4's built-in CLHEP and use of a custom CLHEP:
#set(SYSTEM_CLHEP OFF)

#if( NOT EXISTS ${Geant4_INCLUDE_DIR}/CLHEP ) 
#    set(SYSTEM_CLHEP ON)
#endif()



#ROOT path
# commented out the following line so that this procedure will evolve with ROOT:
# set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${Geant4_INCLUDE_DIR} ${CLHEP_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/src/dss2007 ${CMAKE_CURRENT_BINARY_DIR}/include)





#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/G4SBS*.cc ${PROJECT_SOURCE_DIR}/src/dss2007/*.cc ${PROJECT_SOURCE_DIR}/include/Pythia6_tree.C ${PROJECT_SOURCE_DIR}/include/simc_tree.C)
file(GLOB headers ${CMAKE_CURRENT_BINARY_DIR}/include/gitinfo.hh ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/include/wiser_pion.h ${PROJECT_SOURCE_DIR}/src/cteq/*.h ${PROJECT_SOURCE_DIR}/src/dss2007/*.hh ${PROJECT_SOURCE_DIR}/include/Pythia6_tree.h ${PROJECT_SOURCE_DIR}/include/simc_tree.h)

file(GLOB cteqsources ${PROJECT_SOURCE_DIR}/src/cteq/*.c  )
file(GLOB cteqheaders ${PROJECT_SOURCE_DIR}/src/cteq/*.h )

file(GLOB libsources ${PROJECT_SOURCE_DIR}/src/G4SBSRunData.cc ${PROJECT_SOURCE_DIR}/src/G4SBSTextFile.cc)
file(GLOB libheaders ${CMAKE_CURRENT_BINARY_DIR}/include/gitinfo.hh ${PROJECT_SOURCE_DIR}/include/*.hh )

# Generate the dictionary for G4SBS classes deriving from TOb ect:
ROOT_GENERATE_DICTIONARY(
    G__g4sbsroot                   
    G4SBSRunData.hh         
    G4SBSTextFile.hh
    OPTIONS -I${PROJECT_SOURCE_DIR}/include
    MODULE g4sbsroot
    LINKDEF "include/G4SBSLinkDef.h"  # ROOT linkDef file
)
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(sbscteq STATIC ${cteqsources} ${cteqheaders})
SET_TARGET_PROPERTIES(sbscteq PROPERTIES COMPILE_FLAGS "-DCTEQ_TBL_PATH=\\\"${PROJECT_SOURCE_DIR}/cteq-tbls\\\"")

add_executable(g4sbs g4sbs.cc ${sources} ${headers})
add_library(g4sbsroot SHARED ${libsources} ${libheaders} G__g4sbsroot.cxx)

target_link_libraries(g4sbs g4sbsroot ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} sbscteq )
target_link_libraries(g4sbsroot ${ROOT_LIBRARIES} )

#if( ${SYSTEM_CLHEP} )
# check for absence of GEANT4's built-in CLHEP. If it isn't there,
# we have to find it:
if( NOT ${Geant4_builtin_clhep_FOUND})
    find_package(CLHEP REQUIRED)

    include_directories(${CLHEP_INCLUDE_DIR})
    target_link_libraries(g4sbs CLHEP::CLHEP)
    target_link_libraries(g4sbsroot CLHEP::CLHEP)
endif()

option(WITH_G4SBS_GDML "Build g4sbs with GDML output" OFF)
if( WITH_G4SBS_GDML)
    set_target_properties(g4sbs 
    PROPERTIES COMPILE_FLAGS "-DG4SBS_USE_GDML")
endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build g4sbs. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
# need to rework this so that these scripts get installed into a directory that is pointed to by an environment variable that will be set up by a shell script to be 
# sourced to the user's environment.
# both this and the following section copying the ROOT macros to the build directory are now basically obsolete: 
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
  scripts/genrp_4.5GeV2.mac
  scripts/gep_5GeV2.mac
  scripts/gep_8GeV2.mac
  scripts/gep_12GeV2.mac
  #scripts/sidis_template.mac
  #scripts/make_SIDIS_scripts.C
  scripts/vis_sidis.mac
  scripts/sidis_background.mac
  scripts/vis_gep.mac
  scripts/vis_gmn.mac
  scripts/vis_a1n.mac
  scripts/vis_genrp.mac
  scripts/preinit_ckov_scint.mac
  scripts/preinit_ckov_noscint.mac
  scripts/preinit_nockov_scint.mac
  scripts/preinit_nockov_noscint.mac
  scripts/preinit_ckov_noscint_nocalorimeters.mac
  scripts/sidis_example_script.mac
  scripts/default_thresholds.mac
  database/ecal_map.txt
  scripts/transversity.mac
  scripts/a1n_bigbite30deg_disgen.mac
  scripts/default_thresholds.mac
  scripts/sidis_template.mac
  )

foreach(_script ${G4SBS_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

set(G4SBS_ROOT_MACROS
root_macros/C16_analysis.C
root_macros/C16_tree.C
root_macros/C16_tree_with_flux.C
root_macros/DIS_rates_DIS.C
root_macros/DIS_rates_PYTHIA.C
root_macros/ECAL_background_dose_rate.C
root_macros/Fit_SBS_optics_and_spin.C
root_macros/G4SBS_optics_fit.C
root_macros/G4SBS_spin_transport.C
root_macros/Optics_A1n_SIDIS.C
root_macros/Pythia6_minbias_gen.C
root_macros/Pythia6_tree.C
root_macros/RICH_analysis.C
root_macros/SIDIS_DST.C
root_macros/SIDIS_reduced_data_file.C
root_macros/TriggerLogicECalGEP.C
root_macros/a1n_tree.C
root_macros/background_rates_GMN.C
root_macros/flux.C
root_macros/g4sbs_a1n_tree.C
root_macros/g4sbs_tree.C
root_macros/gep_optics_tree.C
root_macros/gep_tree.C
root_macros/gep_tree_July2015.C
root_macros/gep_tree_elastic.C
root_macros/gep_tree_with_spin.C
root_macros/gep_trigger_analysis.C
root_macros/gep_trigger_analysis_elastic.C
root_macros/gep_trigger_analysis_elastic_L2.C
root_macros/gep_trigger_analysis_pythia.C
root_macros/gep_trigger_analysis_pythia_L2.C
root_macros/gmn_tree.C
root_macros/rate_table.C
root_macros/rate_table_DST.C
root_macros/C16_tree.h
root_macros/C16_tree_with_flux.h
root_macros/Pythia6_tree.h
root_macros/SIDIS_DST.h
root_macros/a1n_tree.h
root_macros/g4sbs_a1n_tree.h
root_macros/g4sbs_tree.h
root_macros/gep_optics_tree.h
root_macros/gep_tree.h
root_macros/gep_tree_July2015.h
root_macros/gep_tree_elastic.h
root_macros/gep_tree_with_spin.h
root_macros/gmn_tree.h
)

foreach(_macro ${G4SBS_ROOT_MACROS} )
	       configure_file( ${PROJECT_SOURCE_DIR}/${_macro}
	       ${PROJECT_BINARY_DIR}/${_macro} COPYONLY )
endforeach()

configure_file( ${PROJECT_SOURCE_DIR}/rootlogon_builddir.C ${PROJECT_BINARY_DIR}/rootlogon.C )
configure_file( ${PROJECT_SOURCE_DIR}/rootlogon_install.C ${PROJECT_BINARY_DIR}/root_macros/rootlogon.C )

set(DSS2007_GRIDS HLO.GRID HNLO.GRID KALO.GRID KNLO.GRID PILO.GRID PROLO.GRID PLO.GRID PNLO.GRID)

foreach(_grid ${DSS2007_GRIDS})
  configure_file( ${PROJECT_SOURCE_DIR}/src/dss2007/${_grid} ${PROJECT_BINARY_DIR}/${_grid} COPYONLY )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(G4SBS DEPENDS g4sbs)

# Update header file
# First, try to find a suitable version of python
# Either version 2 or anything >= 3.5 works
find_package( Python2 )
if(${Python2_FOUND})
  set( GIT_INFO_SCRIPT ${PROJECT_SOURCE_DIR}/pullgitinfo.py)
  set( Python_EXECUTABLE ${Python2_EXECUTABLE})
else() ## Look for python >= 3.5
  find_package( Python ) ## The default is to find 3 and above
  if(${Python_FOUND} AND ${Python_VERSION_MAJOR} GREATER_EQUAL 3 AND
      ${Python_VERSION_MINOR} GREATER_EQUAL 5)
      set( GIT_INFO_SCRIPT ${PROJECT_SOURCE_DIR}/pullgitinfo3.py)
  else() ## No supported version found
    message(FATAL_ERROR "No suitable version of python found. Need either python2 or python >= 3.5, but found version ${Python_VERSION}")
  endif()
endif()

add_custom_target( _gitinfo ALL COMMAND ${Python_EXECUTABLE} ${GIT_INFO_SCRIPT} ${PROJECT_SOURCE_DIR} COMMENT "Creating gitinfo.hh" )

#add_custom_command( TARGET g4sbsroot PRE_BUILD COMMAND ${PROJECT_SOURCE_DIR}/pullgitinfo.py ${PROJECT_SOURCE_DIR} COMMENT "Creating gitinfo.hh")

add_dependencies(g4sbs _gitinfo)
add_dependencies(g4sbsroot _gitinfo)

# copy environment setup script(s) for sh and csh to 
configure_file( ${PROJECT_SOURCE_DIR}/g4sbs.sh ${CMAKE_CURRENT_BINARY_DIR}/g4sbs.sh )
configure_file( ${PROJECT_SOURCE_DIR}/g4sbs.csh ${CMAKE_CURRENT_BINARY_DIR}/g4sbs.csh @ONLY)

file(GLOB dss_grids ${PROJECT_SOURCE_DIR}/src/dss2007/*.GRID)

#file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/run_g4sbs_here)

configure_file(${PROJECT_SOURCE_DIR}/g4sbs.rootrc ${CMAKE_CURRENT_BINARY_DIR}/.rootrc)



#TODO: figure out how to get the ROOT_GENERATE_DICTIONARY command to do this automatically:
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libg4sbsroot_rdict.pcm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libg4sbsroot.rootmap DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/G__g4sbsroot.cxx DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
# 

install(TARGETS g4sbs RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS sbscteq ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(TARGETS g4sbsroot LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/g4sbs.sh DESTINATION ${CMAKE_INSTALL_BINDIR})
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/g4sbs.csh DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY ${PROJECT_SOURCE_DIR}/scripts DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY ${PROJECT_SOURCE_DIR}/root_macros DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY ${PROJECT_SOURCE_DIR}/farmsub_scripts DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY ${PROJECT_SOURCE_DIR}/integration_tests DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY ${PROJECT_SOURCE_DIR}/database DESTINATION ${CMAKE_INSTALL_PREFIX}/run_g4sbs_here)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/transversity_grids DESTINATION ${CMAKE_INSTALL_DATAROOTDIR})
#install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/run_g4sbs_here ${CMAKE_INSTALL_PREFIX})
install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${dss_grids} DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/DSS2007_GRIDS)
install(FILES ${PROJECT_SOURCE_DIR}/vis.mac ${PROJECT_SOURCE_DIR}/run_example.mac DESTINATION ${CMAKE_INSTALL_PREFIX}/scripts)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/.rootrc DESTINATION ${CMAKE_INSTALL_PREFIX}/run_g4sbs_here)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/root_macros/rootlogon.C DESTINATION ${CMAKE_INSTALL_PREFIX}/root_macros)

option(INSTALL_SBS_FIELD_MAPS "Download and/or install field maps" OFF)
if(INSTALL_SBS_FIELD_MAPS)
  # configure download of field maps (BigBite, SBS GEP 12 GeV^2 TOSCA, GMN 13.5 GeV^2 TOSCA):
  message("Downloading BigBite field map...")
  file(DOWNLOAD https://hallaweb.jlab.org/12GeV/SuperBigBite/downloads/map_696A.dat ${CMAKE_CURRENT_BINARY_DIR}/map_696A.dat SHOW_PROGRESS)
  message("Downloading SBS TOSCA map for GEP, 12 GeV^2")
  file(DOWNLOAD https://hallaweb.jlab.org/12GeV/SuperBigBite/downloads/GEP_12map0_newheader.table ${CMAKE_CURRENT_BINARY_DIR}/GEP_12map0_newheader.table SHOW_PROGRESS)
  message("Downloading SBS+BB TOSCA map for GMN, 13.5 GeV^2")
  file(DOWNLOAD https://hallaweb.jlab.org/12GeV/SuperBigBite/downloads/GMN13_field_newheader.table ${CMAKE_CURRENT_BINARY_DIR}/GMN13_field_newheader.table SHOW_PROGRESS)

#  file(DOWNLOAD https://hallaweb.jlab.org/12GeV/SuperBigBite/downloads/GEN-RP-tosca_newheader.table ${CMAKE_CURRENT_BINARY_DIR}/GEN-RP-tosca_newheader.table SHOW_PROGRESS)

   message("Downloading Portable BigBite field map generated from GMN global map")		     
  file(DOWNLOAD https://hallaweb.jlab.org/12GeV/SuperBigBite/downloads/BigBitePortableFieldMap.table ${CMAKE_CURRENT_BINARY_DIR}/BigBitePortableFieldMap.table SHOW_PROGRESS)

  message("Downloading portable SBS field map generated from GMN global map, with front clamp only, no pole shims")
  file(DOWNLOAD https://hallaweb.jlab.org/12GeV/SuperBigBite/downloads/SBSPortableFieldMap_FrontClampOnly_NoPoleShims.table ${CMAKE_CURRENT_BINARY_DIR}/SBSPortableFieldMap_FrontClampOnly_NoPoleShims.table SHOW_PROGRESS)

  message("Downloading portable SBS field map generated from GEP global map, with front and rear clamps and pole shims")
  file(DOWNLOAD https://hallaweb.jlab.org/12GeV/SuperBigBite/downloads/SBSPortableFieldMap_TwoClamps_PoleShims.table ${CMAKE_CURRENT_BINARY_DIR}/SBSPortableFieldMap_TwoClamps_PoleShims.table SHOW_PROGRESS)

  message("Downloading GEN-RP field map from Bogdan, New Year's Day 2024")
  file(DOWNLOAD https://hallaweb.jlab.org/12GeV/SuperBigBite/downloads/GEn-RP-SBS_field.table ${CMAKE_CURRENT_BINARY_DIR}/GEn-RP-SBS_field.table SHOW_PROGRESS)	

  message("Other field maps forthcoming... Anybody know how to get CMake to download files AFTER the configure step, during the build step?")

  file(GLOB sbs_field_maps ${CMAKE_CURRENT_BINARY_DIR}/*.table ${CMAKE_CURRENT_BINARY_DIR}/map_696A.dat)

  install(FILES ${sbs_field_maps} DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/fieldmaps)
endif()
