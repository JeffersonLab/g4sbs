# g4sbs
[![Build Status](https://travis-ci.org/JeffersonLab/g4sbs.svg?branch=master)](https://travis-ci.org/JeffersonLab/g4sbs)

Seamus Riordan  
riordan@jlab.org  
March 11, 2016  

## Requirements:

- Geant4 version 10.1 or later
- cmake >= 3.9
- root 5.34 (root version 6 strongly recommended; ROOT 5 no longer actively supported at JLab)
- python
- git (optional)
## For detailed documentation:
See the documentation wiki in the Hall A homepage:

https://hallaweb.jlab.org/wiki/index.php/Documentation_of_g4sbs#Documentation_of_g4sbs:_Overview

## Building and running G4SBS
Build using the standard `cmake` facilities. CTEQ tables built in automatically for DIS
If you would like to get on a mail notification list, contact Andrew Puckett (andrew.puckett@uconn.edu)

**IF YOU ARE DOING DEVELOPMENT, PLEASE DO NOT WORK IN THE master BRANCH**

You need the BigBite field map for this as well, which can be found

http://hallaweb.jlab.org/12GeV/SuperBigBite/downloads/map_696A.dat

To build immediately after cloning (**not in g4sbs/**), create a local "build" directory parallel to the source directory. You can run g4sbs directly in the build folder. 

```shell
mkdir build
cd build
cmake ../g4sbs
#  see below if you get an error about old versions
make
```
To run:
```shell
./g4sbs run_example.mac
```

run_example.mac demonstrates most of the features by example


**NOTE:  On the JLab farm, a newer version of cmake is available**  
**So, if it complains about being too old, try:**
```shell
/apps/cmake/bin/cmake ../g4sbs  
```
#########################################################

## Coordinate systems in hall:

+z is down the nominal beam axis  
+y is "up" (away from gravity)  
+x makes a right handed coordinate system  

## Coordinate systems for detectors:

+z is nominally down the "central axis" in the particle direction  
+x is "down" (into the floor)  
+y makes a right handed coordinate system  

#########################################################

Hits in the GEM chambers have tracking performed on them doing
straight line fits with resolutions.

At the moment tracking is done in a very dumb way for simplicity!  
It will end up averaging all the hits together if there is more 
than one track.  It would be wise to ensure that tr.nhit <= 4

## ROOT file output structure:

```shell
ev.*  
  count # Counts given the beam time/luminosity  
  rate # Counts per second given luminosity  
  solang # Integrating over this variable will give solid angle   
  sigma # Cross section [cm^-2]  
  W2 # Invariant mass squared [GeV^2]  
  xbj # Bjorken-x  
  Q2 # Q2 [GeV^2]  
  th # Polar angle of electron with zaxis [rad]  
  ph # Azimuthal angle of electron with zaxis [rad]  
  Aperp # perp component of asymmetry  
  Apar  # parallel component of asymmetry  
  vx,y,z # Vertex position [cm]  
  ep # Scattered electron momentum [GeV]  
  np # Scattered nucleon momentum [GeV]  
  epx,y,z Lab components of electron momentum  [GeV]  
  npx,y,z Lab components of nucleon momentum  [GeV]  
  nth  # Polar angle of nucleon wrt z-axis [rad]  
  nph     Azimuthal angle of nucleon wrt z-axis [rad]  
  pmperp # Missing momentum perp-component [GeV]  
  pmpar  # Missing momentum par-component [GeV]  
  nucl # Nucleon when scattering type, 0 for neutron, 1 for proton  
  fnucl # Final nucleon type (pion prod may change flavor), 0 for neutron, 1 for proton  
```
-----------------------------------------------
```shell
tr.*
  x # Track x coordinate intercept with z= 0 plane [m]
  y # Track y coordinate intercept with z= 0 plane [m]
  xp # Track dx/dz
  yp # Track dy/dz
  tx,ty,typ,txp
   # "True" track variables defined by track projection at first chamber
  hcal # 0 if no hit in HCAL,  1 if hit
  bb # 0 if no hit in BB cal, 1 if hit
  gemtr # 0 if no track found in GEMs, 1 if track found
  hcx # HCAL x hit position [cm]
  hcy # HCAL y hit position [cm]
  bcx # BB cal x hit position [cm]
  bcy # BB cal y hit position [cm]
  hct # HCAL time-of-flight [ns]
  hctex # HCAL expected time-of-flight from momentum-transfer [ns]
  hclx,y,z
   # HCAL hit position in lab coordinates [cm]
  hcdang # HCAL angular difference between q and nucleon p [rad]
```
-----------------------------------------------
```shell
gen.*
  thbb # BigBite angle [rad]
  thhcal # HCAL angle [rad]
  dbb # BigBite distance from target [m]
  dhc # HCAL distance from target [m]
  Ebeam # Beam energy [GeV]
```

-----------------------------------------------
```shell
ht.*
  ndata # Number of hits in ht.* array
  gid # GEM ID (counting starts at 1)
  x,y,z # GEM hit position (incl. resolution effects) [m]
  dx,dy # dx/dz and dy/dz of track at hit
  tx,ty # "True" GEM hit position (perfect resolution) [m]
```



