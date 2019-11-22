#!/bin/csh

#script to set up the environment for g4sbs:
setenv G4SBS @CMAKE_INSTALL_PREFIX@

if( ! ${?PATH} ) then
    setenv PATH $G4SBS/bin
else
    setenv PATH $G4SBS/bin:${PATH}
endif

set OS=`uname -s`


if( "$OS" == "Darwin" ) then # Mac OS: set DYLD_LIBRARY_PATH to library directory:
    if(! ${?DYLD_LIBRARY_PATH}) then
	setenv DYLD_LIBRARY_PATH $G4SBS/lib
    else
	setenv DYLD_LIBRARY_PATH $G4SBS/lib:${DYLD_LIBRARY_PATH}
    endif
endif

# set LD_LIBRARY_PATH regardless of OS:
if( ! ${?LD_LIBRARY_PATH}) then
    setenv LD_LIBRARY_PATH $G4SBS/lib
else
    setenv LD_LIBRARY_PATH $G4SBS/lib:${LD_LIBRARY_PATH}
endif


