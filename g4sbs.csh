#!/bin/csh

#script to set up the environment for g4sbs:
setenv G4SBS @CMAKE_INSTALL_PREFIX@

if( ! ${?PATH} ) then
    setenv PATH @CMAKE_INSTALL_FULL_BINDIR@
else
    setenv PATH @CMAKE_INSTALL_FULL_BINDIR@:${PATH}
endif

set OS=`uname -s`


if( "$OS" == "Darwin" ) then # Mac OS: set DYLD_LIBRARY_PATH to library directory:
    if(! ${?DYLD_LIBRARY_PATH}) then
	setenv DYLD_LIBRARY_PATH @CMAKE_INSTALL_FULL_LIBDIR@
    else
	setenv DYLD_LIBRARY_PATH @CMAKE_INSTALL_FULL_LIBDIR@:${DYLD_LIBRARY_PATH}
    endif
endif

# set LD_LIBRARY_PATH regardless of OS:
if( ! ${?LD_LIBRARY_PATH}) then
    setenv LD_LIBRARY_PATH @CMAKE_INSTALL_FULL_LIBDIR@
else
    setenv LD_LIBRARY_PATH @CMAKE_INSTALL_FULL_LIBDIR@:${LD_LIBRARY_PATH}
endif


