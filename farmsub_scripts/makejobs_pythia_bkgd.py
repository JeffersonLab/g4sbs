# This script is used to submit g4sbs jobs on the farm. 
# More specifically, this script submit jobs with pythia generator (i.e. Pythia files)
# Example: to run pythia in the sidis configuration (script sidis_pythia_bkgd.mac), type:
# > python2 makejobs_pythia_bkgd.py sidis_pythia_bkgd.mac 10 
# for 10 jobs (obviously, any number of jobs between 1 the number of pythia files you have is good)
# if one wants to include a preinit script (e.g preinit_ckov_noscint.mac to turn off scintillation), type:
# > python2 makejobs_pythia_bkgd.py sidis_pythia_bkgd.mac 10 preinit_ckov_noscint.mac 
# Obviously, before using this script for your own, change all paths indicated by "/your_/path_/to_/your_/..."

#!/usr/bin/python

import sys
import os
import datetime

if len(sys.argv) != 3 and len(sys.argv) != 4:
    print "Supply macro and number of jobs"
    print " you may supply a preinit macro as a last argument"
    sys.exit()

macro = sys.argv[1][:-4]
cwd = os.getcwd() 

if len(sys.argv) == 4:
    macro_pre = sys.argv[3][:-4]

runfiletxt="""#!/bin/csh
#setenv ROOTSYS /apps/root/5.34.13/root
#setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH':'$ROOTSYS
"""


for i in range(0, int(sys.argv[2])):
    mytime = datetime.datetime.today()
    suffix = macro+'_'+mytime.strftime("%Y%m%d_%H%M%S")+'_'+str(i)
    suffix2 = macro+'_'+mytime.strftime("%Y%m%d_%H")
    suffix3 = str(i)

    print 'Creating job ', suffix

    runfile = open("runjob_"+suffix+".sh", 'w')
    runfile.write(runfiletxt)
    runfile.write("\n")
    runfile.write("ln -s /group/exjpsi/eric/g4work/g4sbs-build-new/map_696A.dat .\n")
    runfile.write("ln -s /group/halla/www/hallaweb/html/12GeV/SuperBigBite/downloads/map_696A.dat .\n")
    runfile.write("ln -s /group/halla/www/hallaweb/html/12GeV/SuperBigBite/downloads/GEP_12map0_newheader.table .\n")
    #runfile.write("ln -s /group/halla/www/hallaweb/html/12GeV/SuperBigBite/downloads/GMN13_field.table .\n")
    # not available yet
    runfile.write("ln -s /your_/path_/to_/your_/pythia_files_/Pythia_file_"+suffix3+".root pythia_file.root\n")
    runfile.write("mkdir /volatile/your_/path_/to_/your_/output_/"+suffix2+"\n")
    if len(sys.argv) == 3:
        runfile.write("time ./g4sbs "+macro+".mac > run.out\n")
    if len(sys.argv) == 4:
        runfile.write("time ./g4sbs "+macro_pre+".mac "+macro+".mac > run.out\n")
    runfile.write("mv pythia_bkgd.root /volatile/your_/path_/to_/your_/output_/"+suffix2+"/pythia_bkgd_"+suffix3+".root\n")
# NB: it is important that the name of the root file you produce (in your .mac file) is the same as the first file name in the previous line. If not, your output file will be lost.
    runfile.write("mv run.out /volatile/your_/path_/to_/your_/output_/run_"+suffix+".out\n")
    runfile.close()
    
    subfile = open("g4sbs_work_"+suffix+".jsub", 'w')
    
    subfile.write("PROJECT: sbs\n")
    subfile.write("TRACK: simulation\n")
    #subfile.write("TRACK: debug\n")
    subfile.write("OS: centos7\n")
    subfile.write("JOBNAME: g4sbs_"+suffix+"\n")
    subfile.write("MEMORY: 1500 MB\n") #can be tuned for more memory if needed
#NB: it is recommended to run tests before generating large statistics to evaluate how much memory is needed
    subfile.write("COMMAND: csh runjob_"+suffix+".sh\n")
    subfile.write("OTHER_FILES:\n")
    subfile.write("/your_/path_/to_/your_/submit_dir_/runjob_"+suffix+".sh\n")
    subfile.write("/your_/path_/to_/your_/g4sbs_build_dir_/g4sbs\n")
    subfile.write("/your_/path_/to_/your_/pythia_files_/Pythia_file_"+suffix3+".root\n")
    subfile.write(cwd+'/'+macro+".mac\n")
    if len(sys.argv) == 4:
        subfile.write(cwd+'/'+macro_pre+".mac\n")
    subfile.close()

    os.system("jsub g4sbs_work_"+suffix+".jsub")












