#!/usr/bin/env python

##############################################################################
## Specify a default configuration
conf = {}

##  Specify the environment type. Possible options are:
##    clone:
##      Gets a new copy of the repository. Ideally, this is what one uses when
##      not in the Travis-CI environment. The path to the repository and branch
##      to use can be specified later. Tests already committed changes.
##    travis:
##      Uses the repository set up by Travis-CI.
##    local:
##      Used to test local changes which have not yet been committed. Usually
##      one runs this before any commit just to ensure that the Travis-CI test
##      has a high probavility of success. This will not attempt to modify
##      the local environment of the user.  Specify g4sbs executable and
##      working directory.
conf['env_type']          = 'clone'

## Address of the git repository to clone (only used if env_type set to clone)
conf['g4sbs_repo']        = 'https://github.com/JeffersonLab/g4sbs.git'

## Branch of to checkout when cloning (only used if env_type set to clone)
conf['g4sbs_branch']      = 'cmu_dev'

## Path to the local repository to use (only used if the env_type is not
## set to clone).
conf['local_repo_path']   = '/home/travis/build/Jeffersonlab/g4sbs'

## Path to the build directory of G4SBS (only used if env_type is not
## set to clone).
conf['local_build_path']  = '/home/travis/build/Jeffersonlab/build_g4sbs'

## Path to the field maps directory
conf['fieldmap_dir']      = '/work/halla/sbs/g4sbs_buildtests/fieldmaps'

## Path to temp directory (used only when env_type is set to clone)
conf['tmpdir']            = '/volatile/halla/sbs/g4sbs_buildtests'

## Path to logs directory
conf['logs_dir']          = '/work/halla/sbs/g4sbs_buildtests/logs'

import os
import sys
import re
import glob
import tempfile ## For temporary file/directory creation
import shutil ## To clean up temporary directories created
import subprocess ## So we can run git commands
import datetime

## Internal configuration variables
conf_clone_new_repo = False
conf_build_g4sbs = False
conf_link_fieldmaps = False
conf_make_temp_dir = False
conf_make_build_dir = True
conf_repo_dir = ''
conf_build_dir = ''

## Status of all tests
all_tests_status = True
all_tests_count  = 0
all_tests_passed = 0

## Output log file handle
run_date = datetime.datetime.today();
log_file = ''

## A temporary directory
temp_dir_name = ''

## List of experiments to run
#experiment_script_list = ( 'GEn', 'GMn', 'GEp' )
experiment_script_list = { 'GEn', 'GMn', 'GEp' }
#experiment_script_list = { 'GEp' }

## Exit the program with given code number returned.
## Removes the temporary directory created
def cleanExit(code):
  ## Remove temp dir if it got created
  if conf_make_temp_dir and temp_dir_name:
    shutil.rmtree(temp_dir_name)

  myOutLog('',True)
  myOutLog('',True)
  myOutLogJust('All tests: %3d / %3d passed' %
      (all_tests_passed,all_tests_count))
  if code == 0:
    myOutResult('passed')
  else:
    myOutResult('failed')
  sys.exit(code)

def myLog(text,newline = False):
  if log_file:
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
  ## Log in logfile (if it is already open)
  myLog(text)
def myOutLogJust(text):
  myOut(text.ljust(70,'.'),False)
  myLog(text.ljust(70,'.'),True)

## Print out the result of a test using colors for the terminal
def myOutResult(status):
  status = status.upper()
  color = '' ## Default is no color
  if status == 'DONE':
    color = '\033[94m'
  elif status == 'PASSED':
    color = '\033[92m'
  elif status == 'FAILED':
    color = '\033[91m'
  status = '[' + status + ']'
  pre = ''
  if len(status) < 10:
    pre = pre.rjust(10-len(status),'.')
  myOut(pre + color + status + '\033[m',True)

## Read in a configuration file
def readConfFile(fileName):
  global conf
  f = open(fileName,'r')
  lineCount = 0
  for line in f:
    lineCount += 1
    ## Remove any spacing at the front
    line = line.lstrip()

    ## Remove comments and new line character
    line = line.partition('#')[0].rstrip()

    ## Look for an equal sign in non-empty lines
    if line:
      ## Split line around the equal sign
      var = line.split('=')
      if len(var) != 2:
        myErr('Error reading in configuration %s in line %d: Equal sign missing.'
          % (fileName, lineCount),True)
        cleanExit(3)
      elif var[0] in conf: ## If the key is known, we read it in
        conf[var[0]] = var[1]
      else:
        myErr('Error reading in configuration %s in line %d: Variable %s is unknown.'
          % (fileName, lineCount, var[0]),True)
        cleanExit(3)

## Implement function to parse git commands as suggested by
## http://stackoverflow.com/questions/15079442/running-git-clone-command-in-python
#def git(*args):
#  return subprocess.check_call(['git'] + list(args))
## Use the same technique for ll other commands
def run(cmd,*args):
  return subprocess.check_call([cmd] + list(args),stdout=log_file,stderr=log_file)
def runGetStatus(cmd,*args):
  try:
    run(cmd,*args)
  except subprocess.CalledProcessError:
    return 1
  except OSError:
    return 2
  ## If all is good, return status code 0
  return 0
def runExit(cmd,*args):
  status = runGetStatus(cmd,*args)
  if status==1:
    myErr('\n\n\nError: %s produced an error while running.' % cmd,True)
    myErr('Please check the output and try again. Stopping tests.',True)
    cleanExit(2)
  if status==2:
    myErr('Please check your environment and try again. Stopping tests.',True)
    myErr(('\n\n\nError: %s not found.' % cmd),True)
    cleanExit(2)

def init():
  global temp_dir_name, log_file
  global conf_clone_new_repo, conf_build_g4sbs
  global conf_link_fieldmaps, conf_make_temp_dir
  global conf_repo_dir, conf_build_dir

  ## Some configuration after specifing the environment
  if conf['env_type'].lower() == 'clone':
    conf_clone_new_repo = True
    conf_build_g4sbs = True
    conf_link_fieldmaps = True
    conf_make_temp_dir = True
    conf_make_build_dir = True
  elif conf['env_type'].lower() == 'travis':
    conf_clone_new_repo = False
    conf_build_g4sbs = True
    conf_link_fieldmaps = True
    conf_make_temp_dir = False
    conf_repo_dir = conf['local_repo_path']
    conf_build_dir = conf['local_build_path']
    conf_make_build_dir = True
  elif conf['env_type'].lower() == 'local':
    conf_clone_new_repo == False
    conf_build_g4sbs = False
    conf_link_fieldmaps = False
    conf_make_temp_dir = False
    conf_repo_dir = conf['local_repo_path']
    conf_build_dir = conf['local_build_path']
    conf_make_build_dir = False
  else:
    myErr("Invalid env_type selected")
    sys.exit(1000)

  ## Open the log_file
  log_file = open(conf['logs_dir'] + '/buildtests_' +
      run_date.strftime('%Y%m%d_%H:%M') + '.log','w')

  ## Make temporary directory?
  if conf_make_temp_dir:
    temp_dir_name = tempfile.mkdtemp(dir=conf['tmpdir'])
    myOutLog('Setting up temporary directory: ' + temp_dir_name, True)
    conf_repo_dir = temp_dir_name + '/g4sbs'
    conf_build_dir = temp_dir_name + '/build'

  ## Make the build directory?
  if conf_make_build_dir:
    os.mkdir(conf_build_dir)

## Clone a fresh copy of the repository
def cloneRepo():
  pwd = os.getcwd() ## Get the current working directory
  myOutLog('\n\nSetting up fresh copy of G4SBS from ' + conf['g4sbs_repo'],True)
  myOutLogJust('Cloning')
  runExit('git','clone',conf['g4sbs_repo'],conf_repo_dir)
  os.chdir(conf_repo_dir)
  ## Fetch all remote branches
  runExit('git','fetch')
  ## Switch to branch to test
  runExit('git','checkout',conf['g4sbs_branch'])
  ## Show last commit
  runExit('git','log','-1')
  os.chdir(pwd) ## Change back to the previous working directory
  myOutResult('done')

## Build G4SBS
def buildG4SBS():
  pwd = os.getcwd() ## Get the current working directory

  ## Change into the build directory
  os.chdir(conf_build_dir)

  myOutLog('\n\n')
  myOutLogJust('Building G4SBS')
  ## Run cmake
  runExit('cmake','-DCMAKE_BUILD_TYPE=RelWithDebInfo',conf_repo_dir)
  ## Run make
  runExit('make')
  myOutResult('done')
  os.chdir(pwd) ## Change back to the previous working directory

## Link in Field Maps
def linkFieldMaps():
  ## Change into the repository directory
  pwd = os.getcwd() ## Get current working directory
  os.chdir(conf_repo_dir)
  myOutLog('\n\nLinking in field maps.',True)
  ## Link in all field maps
  helper_files = glob.glob(conf['fieldmap_dir'] + '/*')
  for filename in helper_files:
    myOutLog('Linking ' + filename + ' -> ' + 
        filename.replace(conf['fieldmap_dir'] +'/',''), True)
    os.symlink(filename,filename.replace(conf['fieldmap_dir']+'/',''))
  os.chdir(pwd) ## Change back into previous working directory

## Run the tests on a script and return the status
def runScript(script,opticalTracking = False):
  preinit = '--pre=scripts/preinit_'
  if opticalTracking:
    preinit += 'ckov_scint.mac'
  else:
    preinit += 'nockov_noscint.mac'
  status = runGetStatus(conf_build_dir+'/g4sbs',preinit,'--post='+script)
  if status==0:
    return True
  elif status==1:
    myErr('\n\n\nThere was an error running G4SBS on script ' + script,True)
    return False
  elif status==2:
    myErr('\n\n\nError: ' + conf_build+dir + '/g4sbs not found.',True)
    return False
  ## Otherwise just return false
  return False

## Process the ROOT file produced by the script
def checkRootfile(exp,script,optical):
  rootScript = 'integration_tests/check_rootfile.C'
  ## Check if the rootScript exists, because apparently ROOT is too silly
  ## to be able to return a non-zero exit code if the macro does not exist :(
  if not os.path.isfile(rootScript):
    myLog(rootScript + ' macro does not exist',True)
    return False
  script = script.replace('scripts/','') + '.root'
  rootfile = script.replace('.mac','')
  scriptParams = '(\"' + rootfile + '\",\"' + exp + '\",' + optical +')'
  status = runGetStatus('root','-b','-q','-x',rootScript + scriptParams)
  if status==0:
    return True
  elif status==1:
    myLog('Processing ROOT file for script ' + script + ' failed.',True)
    myLog('Please check the output for details.',True)
    return False
  elif status == 2:
    myErr('root command not found',True)
    cleanExit(4)

  return False

## Main Function (program starts here)
def main(conf_file):
  global all_tests_status, all_tests_count, all_tests_passed

  ## Ensure we catch the keyboard interruptions
  try:
    ## If a config file was specified, read it in
    if conf_file:
      readConfFile(conf_file)

    ## Initialize the environment
    init()

    ## Do we need to clone a fresh copy of the repository?
    if conf_clone_new_repo:
      cloneRepo()
    else:
      myOutLog('Using existing repository under ' + conf_repo_dir,True)

    ## Build G4SBS?
    if conf_build_g4sbs:
      buildG4SBS()
    else:
      myOutLog('Using already compiled G4SBS binary under ' + conf_build_dir,True)

    ## Do we need to link in field maps?
    if conf_link_fieldmaps:
      linkFieldMaps()

    ## Change into the repository directory and beging processing the
    ## sample experiment's scripts
    os.chdir(conf_repo_dir)
    for exp in experiment_script_list:
      scriptlist = glob.glob('scripts/' + exp.lower() + '_*GeV2.mac')
      myOutLog('\n\nFound %d %s macros. Beginning tests (with no optical photons):'
          % (len(scriptlist),exp), True)
      for script in scriptlist:
        all_tests_count += 1
        myLog('\n\n\n')
        myOutLogJust('  Testing ' + script.replace('scripts/',''))
        status = runScript(script,False)
        if status:
          ## Now check ROOT file
          myLog('    \n\n\nChecking ROOT file',True)
          status = checkRootfile(exp,script,'false')

        ## After all tests check the status again
        myLog('  Test Result of ' + script.replace('scripts/',''))
        if status:
          myLog('  PASSED',True)
          myOutResult('passed')
          all_tests_passed += 1
        else:
          myLog('  FAILED!!',True)
          myOutResult('failed')

    ## Did all tests pass?
    if all_tests_count == all_tests_passed:
      cleanExit(0)
    else:
      cleanExit(404)

  except KeyboardInterrupt:
    myOutLog('\n\nUser pressed the keyboard interrupt. Exiting early.',True)
    cleanExit(505)


## Run main() function if this file is run directly in the command line
## (as opposed to being imported)
if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1])
  else:
    main('')
