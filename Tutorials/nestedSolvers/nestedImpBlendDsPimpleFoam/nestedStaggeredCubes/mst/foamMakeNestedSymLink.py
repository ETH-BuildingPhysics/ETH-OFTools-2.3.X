#!/usr/bin/env python3

#======================================
# name:    foamMakeNestedSymLink.py
# author:  Marcel Vonlanthen
# date:    04.12.2013
# contact: vonlanthen@arch.ethz.ch
# purpose: Create symbolic link for a nested CFD
# version: 1.0
#======================================

#======================================
# import modules
#======================================
import os
import sys
#import shutil
import argparse
import csv

#import numpy as np
#import scipy as sp


#======================================
# parse arguments
#======================================
parser = argparse.ArgumentParser(description='Purpose: change the world, dude! '
                                             'And also create symbolic link for a nested CFD.')
parser.add_argument('-nd',
                    '--nestedDir',
                    dest='ntdDir',
                    type=str,
                    required=True,
                    help='path to the nested case (/my/path/whatever/nestedCase)'        
                    )
parser.add_argument('-ts',
                    '--timeStep',
                    dest='ts',
                    type=float,
                    required=True,
                    help='time step used to create the symlink'
                    )
parser.add_argument('-co',
                    '--cleanOnly',
                    dest="cleanOnly",
                    action="store_true",                    
                    help='only remove existing sym link without creating new ones.')
parser.add_argument('-so',
                    '--serialOnly',
                    dest="serialOnly",
                    action="store_true",                    
                    help='only create serial sym link. Existing only are removed before creating the new ones.')
args = parser.parse_args()


#======================================
# function definitions
#======================================  
def getFoamTs(nb):
    '''
    converte a float or int (nb) into a valid OpenFOAM timestep (26 instead of 26.0) and return
    it as nbstr, a string. nbstr can be used in a path for example.
    '''
    nbstr = str(nb)
    if nbstr[-2:]=='.0':        # if nbstr='26.0'
        return nbstr[:-2]
    elif nbstr[-1]!='0':        #if nbstr='2.8' of '1.897' or '1.008': valid OF time
        return nbstr
    elif nbstr[-1]=='0' and len(nbstr)==1:   #if nbstr='0': valid OF time
        return nbstr
    else:                      #if nbrstr='26.000' or '3.4500'             
        trailingZero = True
        while trailingZero==True:
            nbstr = nbstr[:-1]
            if nbstr[-1]=='0':
                pass
            elif nbstr[-1]=='.':
                nbstr = nbstr[:-1]
                trailingZero = False
                return nbstr
            else:
                return nbstr

def createSerialSymLink(folder):
    src = os.path.join(ntdPath,folder)
    tgt = os.path.join(mstPath,folder,ntdBname)
    if os.path.exists(tgt)==False:
        if os.path.exists(src):
            os.symlink(src, tgt)
        else:
            pass
            
    else:
        pass

def createParallelSymLink(folder,processorList):
    for proc in processorList:
        src = os.path.join(ntdPath,proc,folder)
        tgt = os.path.join(mstPath,proc,folder,ntdBname)
        if os.path.exists(tgt)==False:
            if os.path.exists(src):
                os.symlink(src, tgt)
            else:
                pass
        else:
            pass  
        
def removeSymLink(symLink):
    if os.path.islink(symLink):   
        os.unlink(symLink)  
#======================================
# main
#======================================
runDir = os.getcwd()

#print(os.path.exists(args.mstDir))
#print(os.path.exists(args.ntdDir))

ts = getFoamTs(args.ts)
serialOnly = args.serialOnly

# absolute path with dirname and basename
mstPath = os.path.abspath(runDir)
ntdPath = os.path.abspath(args.ntdDir)

if os.path.exists(ntdPath):
    pass
else:
    sys.exit('Path "'+ntdPath+'" does not exist. Check your path. Exit.')
    
# directory name (path to the OpenFOAM case)
mstDname = os.path.dirname(mstPath)
ntdDname = os.path.dirname(ntdPath)
# base name (case directory only)
mstBname = os.path.basename(mstPath)
ntdBname = os.path.basename(ntdPath)

#print(mstPath)
#print(ntdPath)
#print(mstDname)
#print(ntdDname)
#print(mstBname)
#print(ntdBname)


# find all processor folder in master and nested, and check if the decomposition is the same
mstProcessors = [folder for folder in os.listdir(mstPath) if folder.startswith('processor')]
ntdProcessors = [folder for folder in os.listdir(ntdPath) if folder.startswith('processor')]
if len(mstProcessors)!=len(ntdProcessors):
    sys.exit(mstBname+' and'+ntdBname+' must have the same number of processors. Exit.')
elif len(mstProcessors)==0:
    serialOnly = True
    

# remove existing serial sym links
symLink = os.path.join(mstPath,'constant',ntdBname)
removeSymLink(symLink)

symLink = os.path.join(mstPath,'system',ntdBname)    
removeSymLink(symLink)

symLink = os.path.join(mstPath,ts,ntdBname)      
removeSymLink(symLink)

# remove existing parallel sym links
for proc in mstProcessors:
    symLink = os.path.join(mstPath,proc,'constant',ntdBname)
    removeSymLink(symLink)
    symLink = os.path.join(mstPath,proc,ts,ntdBname)
    removeSymLink(symLink)

# create symlink. If a folder or file with a symilar name exists, skip and report
if args.cleanOnly==False:        
    # create serial sym links in constant, system and ts
    createSerialSymLink('constant')
    createSerialSymLink('system')
    createSerialSymLink(ts)
   
    # create parallel sym links
    if serialOnly==False:
        createParallelSymLink('constant',mstProcessors)
        createParallelSymLink(ts,mstProcessors)
