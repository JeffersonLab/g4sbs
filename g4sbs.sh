#!/bin/sh

#script to set up the environment for g4sbs:
export G4SBS=${CMAKE_INSTALL_PREFIX}

if test "x$PATH" = "x" ; then
    export PATH=${CMAKE_INSTALL_FULL_BINDIR}
else
    export PATH=${CMAKE_INSTALL_FULL_BINDIR}:$PATH
fi

OS=`uname -s`


if [ "$OS" == "Darwin" ]
then # Mac OS: set DYLD_LIBRARY_PATH to library directory:
    if test "x$DYLD_LIBRARY_PATH" = "x"; then
	export DYLD_LIBRARY_PATH=${CMAKE_INSTALL_FULL_LIBDIR}
    else
	export DYLD_LIBRARY_PATH=${CMAKE_INSTALL_FULL_LIBDIR}:$DYLD_LIBRARY_PATH
    fi
fi

# set LD_LIBRARY_PATH regardless of OS:
if test "x$LD_LIBRARY_PATH" = "x"; then
    export LD_LIBRARY_PATH=${CMAKE_INSTALL_FULL_LIBDIR}
else
    export LD_LIBRARY_PATH=${CMAKE_INSTALL_FULL_LIBDIR}:$LD_LIBRARY_PATH
fi


