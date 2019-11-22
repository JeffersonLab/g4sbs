#!/bin/sh

#script to set up the environment for g4sbs:
export G4SBS=${CMAKE_INSTALL_PREFIX}

if test "x$PATH" = "x" ; then
    export PATH=$G4SBS/bin
else
    export PATH=$G4SBS/bin:$PATH
fi

OS=`uname -s`


if [ "$OS" == "Darwin" ]
then # Mac OS: set DYLD_LIBRARY_PATH to library directory:
    if test "x$DYLD_LIBRARY_PATH" = "x"; then
	export DYLD_LIBRARY_PATH=$G4SBS/lib
    else
	export DYLD_LIBRARY_PATH=$G4SBS/lib:$DYLD_LIBRARY_PATH
    fi
fi

# set LD_LIBRARY_PATH regardless of OS:
if test "x$LD_LIBRARY_PATH" = "x"; then
    export LD_LIBRARY_PATH=$G4SBS/lib
else
    export LD_LIBRARY_PATH=$G4SBS/lib:$LD_LIBRARY_PATH
fi


