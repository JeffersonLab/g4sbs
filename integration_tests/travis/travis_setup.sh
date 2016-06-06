#!/bin/bash

## Get the path of the repo
repo_dir=`pwd`

www_base=http://emilio.phys.cmu.edu/~cornejo/public/g4sbs
www_comp=${www_base}/compiled
www_src="${www_base}/source"
www_dat=${www_base}/data

cmake_file=cmake-3.5.2-Linux-x86_64.tar.gz
root_file=root_v5.34.32.Linux-ubuntu12-x86_64-gcc4.6.tar.gz
xercesc_file=xerces-c-3.1.3.tar.gz
geant4_file=geant4.10.01.p01.tar.gz
geant4_share=Geant4-10.1.1
geant4_datalist_file=geant4_data.txt
fieldmaps_file=fieldmaps.tar.gz


###############################################################################
## Start the setup
mkdir logs

APPS_DIR=${HOME}/apps

## Get the environment file
if [ ! -e env.sh ]; then
  wget -c ${www_base}/env.sh
fi

## Make the directory structure
cd ${HOME}
mkdir -p ${APPS_DIR}
mkdir g4sbs_buildtests

## Define the installer function

function local_install()
{
  echo "Installing $1 from archive $2"
  cd ${APPS_DIR}
  if [ ! -d $1 ]; then
    wget -c ${www_comp}/$2
    mkdir $1
    cd $1
    tar xf ../$2
    cd ..
    rm $2
  fi
}

local_install cmake ${cmake_file}
local_install xercesc ${xercesc_file}
local_install geant4 ${geant4_file}
local_install root ${root_file}

## Install the Geant4 data files
geant4_version=`basename ${geant4_file} .tar.gz`
cd ${HOME}/apps/geant4/${geant4_version}/share/${geant4_share}
mkdir data
cd data
wget -c ${www_dat}/${geant4_datalist_file}
#geant4_datafiles=(G4NDL.4.5 G4EMLOW.6.48 G4PhotonEvaporation.3.2 G4RadioactiveDecay.4.3.1 G4SAIDDATA.1.1 G4NEUTRONXS.1.4 G4ABLA.3.0 G4PII.1.3 RealSurface.1.0 G4ENSDFSTATE.1.2.1 G4TENDL.1.0)
for dat in `cat ${geant4_datalist_file}`; do
  if [ ! -e ${dat}.installed ]; then
    wget -c ${www_dat}/${dat}.tar.gz
    tar xf ${dat}.tar.gz
    rm ${dat}.tar.gz
    echo "installed" > ${dat}.installed
  fi
  touch ${dat}.installed
done

## Install the g4sbs field maps
cd ${HOME}
if [ ! -d fieldmaps ]; then
  wget -c ${www_base}/${fieldmaps_file}
  tar xvf ${fieldmaps_file}
  rm ${fieldmaps_file}
  ls
  pwd
fi
