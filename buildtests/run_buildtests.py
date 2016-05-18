#!/usr/bin/env python

## Specify directory name with field maps
fieldmap_dir='/work/halla/sbs/g4sbs_buildtests/fieldmaps/'

## Specify address of git repository
g4sbs_repo='https://github.com/JeffersonLab/g4sbs.git'

## Specify path to temp directory
g4sbs_tmpdir='/volatile/halla/sbs/g4sbs_buildtests/'

import os
import sys
import re
import glob
import tempfile ## For temporary file/directory creation
import shutil ## To clean up temporary directories created
import subprocess ## So we can run git commands
import datetime

## Output log file handle
run_date = datetime.datetime.today();
log_file = open('logs/buildtests_' + run_date. strftime('%Y%m%d_%H:%M') + '.log','w')

## A temporary directory
temp_dir_name = tempfile.mkdtemp(dir=g4sbs_tmpdir)

## Exit the program with given code number returned.
## Removes the temporary directory created
def cleanExit(code):
  ## Remove temp dir
  shutil.rmtree(temp_dir_name)
  sys.exit(code)

def myLog(text,newline = False):
  log_file.write(text)
  if newline:
    log_file.write('\n')
  log_file.flush()
def myOut(text,newline = False):
  sys.stdout.write(text)
  if newline:
    sys.stdout.write('\n')
  sys.stdout.flush()
def myOutLog(text,newline = False):
  ## Print to stdout
  myOut(text,newline)
  ## Also log
  myLog(text,newline)
def myErr(text,newline = False):
  ## Print to stderr
  sys.stderr.write(text)
  if newline:
    sys.stderr.write('\n')
  sys.stderr.flush()
  ## Log in logfile
  myLog(text)

## Implement function to parse git commands as suggested by
## http://stackoverflow.com/questions/15079442/running-git-clone-command-in-python
#def git(*args):
#  return subprocess.check_call(['git'] + list(args))
## Use the same technique for ll other commands
def run(cmd,*args):
  return subprocess.check_call([cmd] + list(args),stdout=log_file,stderr=log_file)
def runExit(cmd,*args):
  try:
    run(cmd,*args)
  except subprocess.CalledProcessError:
    myErr('\n\n\nThere was an error running the command %s.' % cmd,True)
    myErr('Please check the output and try again. Stopping tests.',True)
    cleanExit(2)
  except OSError:
    myErr(('\n\n\nError: %s not found.' % cmd),True)
    myErr('Please check your environment and try again. Stopping tests.',True)
    cleanExit(2)

## Main Function (program starts here)
def main():
  os.chdir(temp_dir_name)
  myOutLog('Setting up temporary directory: ' + temp_dir_name)

  ## Clone a fresh copy of the repository
  myOutLog('\n\nCloning fresh copy of G4SBS from ' + g4sbs_repo,True)
  runExit('git','clone',g4sbs_repo)

  ## Build G4SBS
  myOutLog('\n\nBuilding G4SBS.....')
  os.mkdir('build')
  os.chdir('build')
  ## Run cmake
  runExit('cmake','-DCMAKE_BUILD_TYPE=RelWithDebInfo','../g4sbs')
  ## Run make
  runExit('make')
  myOut('Success',True)

  ## Now change back to the G4SBS directory
  myOutLog('\n\nLinking in field maps.',True)
  ## Link in all field maps
  helper_files = glob.glob(fieldmap_dir + '*')
  for filename in helper_files:
    myOutLog('Linking ' + filename + ' -> ' + 
        filename.replace(fieldmap_dir,''), True)
    os.symlink(filename,filename.replace(fieldmap_dir,''))

  ## Generate a pre-init script that will disable optical photons
  f_preinit = open('buildtests_preinit.mac','w')
  f_preinit.write('## Disable optical photons\n')
  f_preinit.write('/g4sbs/usescint false\n')
  f_preinit.write('/g4sbs/useckov false\n')
  f_preinit.close()

  ## Get a list of scripts to run
  scriptlist = glob.glob('../g4sbs/scripts/*GeV2.mac')
  myOutLog('\n\nFound %d macros. Beginning tests (with no optical photons)...'
      % len(scriptlist), True)
  for script in scriptlist:
    myOutLog('  Testing ' + script.replace('../g4sbs/',''))
    runExit('time','./g4sbs','--pre=buildtests_preinit.mac',
        '--post=' + script)
    myOut('.....Success',True)

  ## So far so good, all tests passed.
  myOutLog('Build test successful....',True)

  ## Exit gracefully
  cleanExit(0)

## Run main() function if this file is run directly in the command line
## (as opposed to being imported)
if __name__ == '__main__':
  main()
