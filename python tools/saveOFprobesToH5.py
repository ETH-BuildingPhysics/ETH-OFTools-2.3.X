#! /usr/bin/env python

# -*- coding: utf-8 -*-

#=============================================================================#
# Import modules
#=============================================================================#

import re
import os
import sys
import time
import h5py
import argparse

# math modules
import numpy as np
import scipy as sp
import scipy.linalg as spla
import scipy.signal as spsi
import math as math

# mpl.use('QT4Agg')         #load qt4 backend (optional)

import pyFlowStat.PointProbe as pp
import pyFlowStat.PointProbeFunctions as ppf
import pyFlowStat.TurbulenceTools as tt
import pyFlowStat.TriSurface as ts

#=============================================================================#
# Functions
#=============================================================================#
def is_number(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def loadOFProbeList(relpath):
    absPath = os.path.abspath(relpath)
    if os.path.exists(absPath):
        print('    loading file "'+str(relpath)+'"')
        return pp.getOFPointProbeList(absPath,reshape=True,createDict=False)
    else:
        sys.exit('Path exist: False. Script abort.')


def appendPPList(runList,createDict):
    runs = []
    
    # load all runs
    for run in runList:
        runs.append(loadOFProbeList(run))
           
    # append runs to the first one of the list
    for runi in range(1,len(runs)):
        print('    appending data from file '+str(runi)+' to data from file 0')
        for pti in range(len(runs[runi])):
            runs[0][pti].appendProbe(runs[runi][pti],rmOverlap='self',createDict=False)    
    
    if createDict==True:
        print('    create data dict')
        for pt in runs[0]:        
            pt.createDataDict(action=createDict)
            
    return runs[0]

def getSubDirList(dir):
    return [name for name in os.listdir(dir) if os.path.isdir(os.path.join(dir, name))]


#=============================================================================#
# parse arguments
#=============================================================================#
parser = argparse.ArgumentParser(description=('Convert an OpenFOAM probe file in a hdf5 probe file. ',
                                              'The OpenFOAM probes must be located in the OF default ',
                                              'path (ex: OFcase/postProcessing/probeName/timestep/U)'))
parser.add_argument('-case',
                    dest='case',
                    type=str,
                    required=False,
                    default=os.getcwd(),
                    help='Path of an OpenFOAM case. Default: execution path'
                    )
parser.add_argument('-probename',
                    dest='probename',
                    type=str,
                    required=True,
                    help='name of the probe to convert. Must be located in postProcessing'
                    )
parser.add_argument('-region',
                    dest='region',
                    type=str,
                    required=False,
                    default=str(),
                    help='specifiy a different region name. Default: region0 (the default OpenFOAM region)'
                    )
parser.add_argument('-keyrange',
                    dest='keyrange',
                    type=str,
                    required=False,
                    choices=['raw','full'],
                    default='raw',
                    help=('Specify which datas must be saved in the hdf5. ',
                          'Choices: raw or full.  Default: raw. See ',
                          'pyFlowStat.PointProbeFunctions.savePPlist_hdf5 for more details.')
                    )
parser.add_argument('-var',
                    dest='var',
                    type=str,
                    required=True,
                    nargs='+',
                    help='Specify variables to save in hdf5.'
                    )
parser.add_argument('-rename',
                    dest='rename',
                    type=str,
                    required=False,
                    default='',
                    help=('Specify an extended name for all the h5 files generated. ',
                          'Example -rename foobar --> U_foobar.h5')
                    )
parser.add_argument('-overwrite',
                    dest="overwrite",
                    action="store_true",
                    default=False,
                    help=('if the hdf5 file already exist, delete it and recreate it. ',
                          'Default behavior: the file already exists, nothing is ',
                          'saved and the execution continues with the other variable')
                    )

args = parser.parse_args()


#=============================================================================#
# Main
#=============================================================================#
#check if case exist and get absolute path
case = args.case
caseAbsPath = str()
if os.path.exists(case):
    caseAbsPath = os.path.abspath(case)
else:
    sys.exit('Case "'+case+'" does not exist. Check your input. Exit')


# check if probename exist and get absolute path
probename = args.probename
probenameAbsPath = os.path.join(caseAbsPath,'postProcessing',probename,args.region)
if os.path.exists(probenameAbsPath):
    pass
else:
    sys.exit('Probename "'+probename+'" does not exist. Check your input. Exit')


# get list of time step in probename and list of time step absolute path
dirList = getSubDirList(probenameAbsPath)
tsList = []
tsListNum=[]
for d in dirList:
    if is_number(d)==True:
        tsList.append(d)
        tsListNum.append(float(d))
tsListNum=np.array(tsListNum)
tsList_idx=np.argsort(tsListNum)
tsList=[tsList[idx] for idx in tsList_idx]
tsList_absPath = [os.path.join(probenameAbsPath,tstep) for tstep in tsList]


# get list of all variable in all time steps
varList_raw = []
for i in range(len(tsList)):
    names = [name for name in os.listdir(tsList_absPath[i]) if os.path.isfile(os.path.join(tsList_absPath[i], name))]
    for name in names:
        varList_raw.append(name)
varList = list(set(varList_raw))


# create a dict with the vars of varList has key, and the list of location containing the vars.
varLocDict = dict()
for var in varList:
    varLocDict[var] = []

for i in range(len(tsList)):
    names = [name for name in os.listdir(tsList_absPath[i]) if os.path.isfile(os.path.join(tsList_absPath[i], name))]
    for name in names:
        varLocDict[name].append(os.path.join(tsList_absPath[i],name))

#-----------------------------------------------------------------------------#
# save all probed varaible in a hdf5 file
#-----------------------------------------------------------------------------#
varToSave = []
#if args.var==None:
#    for key in varLocDict.keys():
#        varToSave.append(key)
#else:
#    for key in args.var:
#        if varLocDict.has_key(key)==True:
#            varToSave.append(key)
#        else:
#            sys.exit('variable "'+key+'" does not exist. Check your input. Exit')

for key in args.var:
    if varLocDict.has_key(key)==True:
        varToSave.append(key)
    else:
        sys.exit('variable "'+key+'" does not exist. Check your input. Exit')

print('case:')
print('    '+caseAbsPath)
print('')

print('probe name:')
print('    '+probenameAbsPath)
print('')

print('keyrange:')
print('    '+args.keyrange)
print('')

print('list of variables converted from OpenFOAM to hdf5:')
for var in varToSave:
    print('    '+var)
print('')

for key in varToSave:
    
    # load data in memory
    print('start converting and saving variable '+key+':')
    createDict = bool()
    if args.keyrange=='full':
        createDict = True
    else:
        createDict = False
    keyData = appendPPList(varLocDict[key],createDict=createDict)
    
    # save data in HDF5
    h5TgtFile = str()
    if args.rename=='':
        h5TgtFile = os.path.join(probenameAbsPath, key+'.h5')
    else:
        h5TgtFile = os.path.join(probenameAbsPath, key+'_'+args.rename+'.h5')
            
    
    if os.path.exists(h5TgtFile)==True and args.overwrite==True:
        os.remove(h5TgtFile)
        print('    Saving variable "'+key+'" in "'+h5TgtFile+'"')
        ppf.savePPlist_hdf5(ppList=keyData, hdf5file=h5TgtFile, keyrange=args.keyrange)

    elif os.path.exists(h5TgtFile)==False:
        print('    Saving variable "'+key+'" in "'+h5TgtFile+'"')
        ppf.savePPlist_hdf5(ppList=keyData, hdf5file=h5TgtFile, keyrange=args.keyrange)

    elif  os.path.exists(h5TgtFile)==True and args.overwrite==False:
        print('file"'+h5TgtFile+'" already exists. Execution continues with the next variable.')
        
    print('')

print('')

