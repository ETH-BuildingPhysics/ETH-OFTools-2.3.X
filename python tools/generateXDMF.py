#! /usr/bin/env python

# -*- coding: utf-8 -*-

# * * * * * * * * * * * * * * * Import modules  * * * * * * * * * * * * * * * #

import re
import os
import sys
import h5py
import argparse
import ntpath

# math modules
import numpy as np
import scipy as sp

# * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * * * #
def is_number(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def pathLeaf(path):
    '''
    get trailing file from a path (works on windows and linux)
    '''
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def sortNumStrList(numStrList):
    '''
    Sort a list of number stored as a list of string. StrNumList is the list of
    number as string.
    '''
    #tsListNum=[]
    numFltList = []
    for nb in numStrList:
        if is_number(nb)==True:
            numFltList.append(float(nb))
    numFltList = np.array(numFltList)
    numStrList_idx = np.argsort(numFltList)
    numStrList_sort = [numStrList[idx] for idx in numStrList_idx]
    return numStrList_sort
    
    

def getH5metaData(hdf5Parser,varList):
    '''
    '''
    # create the returned variable
    varsPerTs = dict()
    varDim = dict()
    nFaces = int()
    nPoints = int()
    
    
    # get all the base keys and timestep keys
    allkeys = hdf5Parser.keys()
    tskeys = allkeys
    try:
        tskeys.pop(tskeys.index('mesh'))
    except:
        raise ValueError('no "mesh" entry in the HDF5 file. Exit.')
    tskeys = sortNumStrList(tskeys)
    
    # fill nFaces and nPoints
    nFaces = hdf5Parser['mesh']['faces'].shape[0]
    nPoints = hdf5Parser['mesh']['points'].shape[0]
    
    # fill varPerTs and varDim
    for ts in tskeys:
        varlist = hdf5Parser[ts].keys()
        varlist.pop(varlist.index('time'))
        if len(args.varList)>0:
            varsPerTs[ts] = args.varList
        else:
            varsPerTs[ts] = varlist
        for var in varlist:
            if varDim.has_key(var)==True:
                pass
            else:
                varShape = hdf5Parser[ts][var].shape
                if len(varShape)>1:
                    varDim[var] = hdf5Parser[ts][var].shape[1]
                else:
                    varDim[var] = 1
    
    return varsPerTs, varDim, nFaces, nPoints
 

def writeTopology(fileName,nFaces,nTab=0,baseTab=2):
    '''
    '''
    btab = baseTab*' '
    tab = nTab*btab   
    lines = (''+tab+'<Topology Type="Triangle" NumberOfElements="'+str(nFaces)+'">\n'
             ''+tab+btab+'<DataStructure Dimensions="'+str(nFaces)+' 3" NumberType="Int" Format="HDF">\n'
             ''+tab+2*btab+''+fileName+':/mesh/faces\n'
             ''+tab+btab+'</DataStructure>\n'
             ''+tab+'</Topology>\n')       
    return lines


def writeGeometry(fileName,nPoints,nTab=0,baseTab=2):
    '''
    '''
    btab = baseTab*' '
    tab = nTab*btab 
    lines = (''+tab+'<Geometry GeometryType="XYZ">\n'
             ''+tab+btab+'<DataStructure Dimensions="'+str(nPoints)+' 3" NumberType="Float" Presicion="4" Format="HDF">\n'
             ''+tab+2*btab+''+fileName+':/mesh/points\n'
             ''+tab+btab+'</DataStructure>\n'
             ''+tab+'</Geometry>\n')
    return lines


def writeAttribute(fileName,nPoints,time,varName,varDim,nTab=0,baseTab=2):
    '''
    '''
    varType = str()
    if varDim==1:
        varType = 'Scalar'
    elif varDim==3:
        varType = 'Vector'
    elif varDim==6:
        varType = 'Tensor6'
    elif varDim==9:
        varType = 'Tensor'
    btab = baseTab*' '
    tab = nTab*btab
    lines = (''+tab+'<Attribute Name="'+varName+'" Center="Node" AttributeType="'+varType+'">\n'
             ''+tab+btab+'<DataStructure Format="HDF" DataType="Float" Precision="4" Dimensions="'+str(nPoints)+' '+str(varDim)+'">\n'
             ''+tab+2*btab+''+fileName+':/'+time+'/'+varName+'\n'
             ''+tab+btab+'</DataStructure>\n'
             ''+tab+'</Attribute>\n')
    return lines
   
   
# * * * * * * * * * * * * * Parse arguments * * * * * * * * * * * * * * * * * #
toolDescription = ('Create the xdmf file from a HDF5 file of surfaces. '
                   'Works only with TruSurface HDF5. the HDF5 can be created '
                   'with the tool "saveOFsurfacesToH5.py".')
                   
parser = argparse.ArgumentParser(description=toolDescription)

parser.add_argument('-h5',
                    '--h5file',
                    dest='h5File',
                    type=str,
                    required=True,
                    help='Path to the HDF5 file.'
                    )
                    
parser.add_argument('-vl',
                    '--varlist',
                    dest='varList',
                    type=str,
                    required=False,
                    nargs='+',
                    default=list(),
                    help='Specify variables to include in the xdmf file.'
                    )

args = parser.parse_args()

# * * * * * * * * * * * * * * * * Main * * * * * * * * * * * * * * * * * * * #
workingDir = os.getcwd()

absH5File = os.path.abspath(args.h5File)
h5FilePath = os.path.dirname(absH5File)
h5FileName = pathLeaf(absH5File)

xdmfFileName = h5FileName[:-3]+'.xdmf'
xdmfFileName = os.path.splitext(h5FileName)[0]+'.xdmf'
absXdmfFile = os.path.join(h5FilePath,xdmfFileName)


# create HDF5 parser
fr = h5py.File(absH5File, 'r')
varsPerTs, varDim, nFaces, nPoints = getH5metaData(hdf5Parser=fr,varList=args.varList)
fr.close()

# tabulation for the entire file
baseTab = 2

# tab size for the header and the bottom lines
btab = baseTab*' '

# start writing the file 
fo = open(absXdmfFile, 'w')

# write the header lines
fo.write('<Xdmf>\n')
fo.write(''+btab+'<Domain>\n')
fo.write(''+2*btab+'<Grid Name="FieldData" GridType="Collection" CollectionType="Temporal">\n')
fo.write('\n')

# write the body
allTs = varsPerTs.keys()
allTs = sortNumStrList(allTs)
for tsStr in allTs:
    tab = 2*btab
    fo.write(''+tab+'<Grid GridType="Collection" CollectionType="Spatial">\n')
    fo.write(''+tab+'<Time Type="Single" Value="'+tsStr+'" />\n')
    fo.write(''+btab+tab+'<Grid Name="'+tsStr+'" Type="Uniform">\n')    
    fo.writelines(writeTopology(fileName=h5FileName,nFaces=nFaces,nTab=4,baseTab=baseTab)) 
    fo.writelines(writeGeometry(fileName=h5FileName,nPoints=nPoints,nTab=4,baseTab=baseTab))
    for var in varsPerTs[tsStr]:
        fo.writelines(writeAttribute(fileName=h5FileName,
                                     nPoints=nPoints,
                                     time=tsStr,
                                     varName=var,
                                     varDim=varDim[var],
                                     nTab=4,
                                     baseTab=baseTab))
    fo.write(''+btab+tab+'</Grid>\n')
    fo.write(''+tab+'</Grid>\n') 
    fo.write('\n')

# write the bottom lines
fo.write('\n')
fo.write(''+2*btab+'</Grid>\n')
fo.write(''+btab+'</Domain>\n')
fo.write('</Xdmf>\n')

fo.close()

# * * * * * * * * * * * * * * * * End * * * * * * * * * * * * * * * * * * * * #