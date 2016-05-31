#!/usr/bin/python

## import the system module
import sys
import re

## Configure the script
outputPath = "../"

## List of experiments to produce sample scripts
exps = [ 'GMn', 'GEn', 'GEp' ]

## Define standard numberOfEvents
numOfEvents = '100'

## Filehandle that will hold kinematics
fOut = 0

## Dictionary that contains kinematic point parameters (and comments)
status = 1
kin = {}
com = {}

## A variable that will hold the column information
cols = [[]]

## Clear a dicionary variable
def clearDictionaryVar(cmd):
  kin[cmd] = 'Undefined'
  com[cmd] = 'Undefined'

## Clears all variables
def clearAll():
  global status,cols,exps,numOfEvents,kin,com,cols
  status = 1
  cols = [[]]
  kin['run']          =  numOfEvents
  clearDictionaryVar('q2')
  clearDictionaryVar('beamcur')
  clearDictionaryVar('target')
  clearDictionaryVar('kine')
  clearDictionaryVar('targpres')
  clearDictionaryVar('targlen')
  clearDictionaryVar('gemres')
  clearDictionaryVar('rasterx')
  clearDictionaryVar('rastery')
  clearDictionaryVar('cerdist')
  clearDictionaryVar('cerdepth')
  clearDictionaryVar('gemsep')
  clearDictionaryVar('bbcaldist')
  clearDictionaryVar('gemconfig')
  clearDictionaryVar('bbfield')
  clearDictionaryVar('48d48field')
  clearDictionaryVar('runtime')
  clearDictionaryVar('beamE')
  clearDictionaryVar('bbang')
  clearDictionaryVar('bbdist')
  clearDictionaryVar('hcaldist')
  clearDictionaryVar('48D48dist')
  clearDictionaryVar('sbsang')
  clearDictionaryVar('thmin')
  clearDictionaryVar('thmax')
  clearDictionaryVar('phmin')
  clearDictionaryVar('phmax')
  clearDictionaryVar('sbsclampopt')
  clearDictionaryVar('exp')
  clearDictionaryVar('sbsmagfield')
  clearDictionaryVar('filename')

def g4sbsCmdPrint(cmd,comment = ''):
  global status,cols,exps,numOfEvents,kin,com,cols
  if cmd in kin and kin[cmd] != 'Undefined':
    if cmd in com and com[cmd] != 'Undefined':
      if comment == '':
        comment = "## " + com[cmd];
      else:
        comment = "## " + com[cmd] + "(" + comment + ")"
    elif comment != '':
        comment = "## " + comment
    if comment != '':
      fOut.write("/g4sbs/%s %s   %s\n" % (cmd.ljust(15),kin[cmd].ljust(15),comment))
    else:
      fOut.write("/g4sbs/%s %s\n" % (cmd.ljust(15),kin[cmd]))


def makeKin():
  global status,cols,exps,numOfEvents,kin,com,cols,fOut
  print("%s" % kin['q2']),
  fOut = open('%s%s_%sGeV2.mac' % (outputPath,kin['exp'].lower(),kin['q2']),'w')
  fOut.write('## Configure G4SBS for %s (Q^2 = %s GeV^2)\n' % (kin['exp'],kin['q2']))
  kin['filename'] = "%s_%sGeV2.root" % (kin['exp'].lower(),kin['q2']);
  g4sbsCmdPrint('filename',"Output rootfile")

  fOut.write('\n## Configure Experiment\n')
  g4sbsCmdPrint('exp')

  fOut.write("\n## Configure the target\n")
  g4sbsCmdPrint('target')
  g4sbsCmdPrint('targpres','Target pressure')
  g4sbsCmdPrint('targlen','Target Length')

  fOut.write('\n## Configure generator settings\n')
  g4sbsCmdPrint('kine','Generator')
  g4sbsCmdPrint('runtime')
  g4sbsCmdPrint('beamcur')
  g4sbsCmdPrint('rasterx')
  g4sbsCmdPrint('rastery')
  g4sbsCmdPrint('beamE')
  g4sbsCmdPrint('thmin')
  g4sbsCmdPrint('thmax')
  g4sbsCmdPrint('phmin')
  g4sbsCmdPrint('phmax')

  fOut.write('\n## Configure standard detector settings\n')
  g4sbsCmdPrint('gemres')
  g4sbsCmdPrint('cerdist')
  g4sbsCmdPrint('cerdepth')
  g4sbsCmdPrint('gemsep')
  g4sbsCmdPrint('bbcaldist')
  g4sbsCmdPrint('gemconfig')
  g4sbsCmdPrint('hcaldist')
  g4sbsCmdPrint('hcalvoffset')
  g4sbsCmdPrint('sbsclampopt')

  fOut.write('\n## Configure the magnets\n')
  g4sbsCmdPrint('bbfield')
  g4sbsCmdPrint('48d48field')
  g4sbsCmdPrint('bbang')
  g4sbsCmdPrint('bbdist')
  g4sbsCmdPrint('sbsang')
  g4sbsCmdPrint('48D48dist')
  g4sbsCmdPrint('sbsmagfield')

  fOut.write('\n## Run %s events\n' % numOfEvents)
  g4sbsCmdPrint('run')
  fOut.close()

## Define the main function
def main():
  global status,cols,exps,numOfEvents,kin,com,cols
  for exp in exps:
    clearAll()
    kin['exp'] = exp.lower()
    print("Making sample scripts for %s for Q^2=(" % exp),
    f = open('%s_kinematics.txt' % exp.lower(), 'r')
    for line in f:
      ## Remove comments and newline character
      line = line.partition('#')[0].rstrip()

      ## Look for a configuration line
      if re.search(r'^config\.',line):
        ## Remove config. from the line, and split by equal sign to get
        ## the variable name, and the setting
        var = line[7:].split('=')
        kin[var[0]] = var[1]

      ## Look for a configuration comment line
      elif re.search(r'^config_comment\.',line):
        var = line[15:].split('=')
        com[var[0]] = var[1]

      ## Look for Table configuration
      elif re.search(r'^config_table=',line):
        var = line[13:].split(',')
        for i in range(0,len(var)):
          col = var[i].split(':')
          if len(col) == 1:
            col.append('')

          ## Expand column matrix
          if i == 0:
            cols[0].extend([col[0],col[1]])
          else:
            cols.append([col[0],col[1]])
        del col

      ## Now for the remaining non-empty lines, parse them in accordance to the
      ## table configuration that was read in a previous line
      elif line:
        ## Split them up by a space and strip all spaces
        row = line.split()
        if len(row) != len(cols):
          print 'Number of entries in this line does not match number of entries in table configuration.'
          status = 0
        else:
          for i in range(len(row)):
            ## If the column is bdl then convert to magnetic field in Tesla
            if cols[i][0] == 'bdl':
              row[i] = str(round(float(row[i])/1.2192,2))
              cols[i][0] = 'sbsmagfield'
              cols[i][1] = 'tesla'
            if cols[i][1] != '':
              kin[cols[i][0]] = row[i] + ' ' + cols[i][1]
            else:
              kin[cols[i][0]] = row[i]
        ## Print/write out kineamtic configuration if all is well
        if status:
          makeKin()
    f.close()
    print(") GeV^2")

## Run main() function if this file is run directly in the command line
## (as opposed to being imported)
if __name__ == '__main__':
  main()
