APPS_DIR=${HOME}/apps

## cmake
export CMAKE_VERSION=3.5.2
export CMAKE_HOME=${APPS_DIR}/cmake/cmake-${CMAKE_VERSION}-Linux-x86_64
export PATH=${CMAKE_HOME}/bin:${PATH}

## XercesC
export XERCESC_VERSION=3.1.3
export XERCESC_HOME=${APPS_DIR}/xercesc/${XERCESC_VERSION}
export XERCESC_INCLUDE_DIR=${XERCESC_HOME}/include
export LD_LIBRARY_PATH=${XERCESC_HOME}/lib:${LD_LIBRARY_PATH}
export PATH=${XERCESC_HOME}/bin:${PATH}

## geant4
export GEANT4_VERSION=4.10.01.p01
source ${APPS_DIR}/geant4/geant${GEANT4_VERSION}/bin/geant4.sh

## Setup ROOT
export ROOT_VERSION=v5.34.32
source ${APPS_DIR}/root/${ROOT_VERSION}/bin/thisroot.sh
