/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application

Description

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "list"
#include "ListListOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    labelList allIndices;
    labelListList procIndices;
    allIndices.setSize(14,0);

    label indicesPerProc = 0;

    procIndices.setSize(Pstream::nProcs());

    //procIndices[Pstream::myProcNo()]=Pstream::myProcNo();

    if(Pstream::master())
    {
        forAll(allIndices,I)
        {
            allIndices[I]=I;
        }

        indicesPerProc = floor(scalar(allIndices.size())/(Pstream::nProcs()-1));
        Info << indicesPerProc << " rest: " << allIndices.size()-(indicesPerProc*(Pstream::nProcs()-1)) << endl;

    }

    Pstream::scatter(allIndices);
    Pstream::scatter(indicesPerProc);

    label start=Pstream::myProcNo()*indicesPerProc;
    label size=indicesPerProc;
    if ((start+size)>allIndices.size())
    {
        size=allIndices.size()-start;
    }
    Pout << "start: " << start << " size: " << size << endl;

    procIndices[Pstream::myProcNo()]=SubList<label>(allIndices,size,start);

    forAll(procIndices[Pstream::myProcNo()],I)
    {
            procIndices[Pstream::myProcNo()][I]=Pstream::myProcNo();
    }

    Pstream::gatherList(procIndices);
    allIndices = ListListOps::combine<labelList>(procIndices,accessOp<labelList>());


    Info << "procIndices " << procIndices << endl;
    Info << "allIndices " << allIndices << "size " << allIndices.size() << endl;
    return 0;

    /*
     * The function I was looking for was:
    listCombineGather() and listCombineScatter().

    But this only works for list, not for arrays... but I guess is not hard to change arrays to lists. Just change the definition is enough. List has index operations.

    For the reference of those who used to program in raw MPI styles, the scatter(), gather(), combineScatter(), and combineGather() are different from the MPI_scatter and MPI_gather, but corresponding to MPI_Bcast and MPI_reduce() respectively.

    MPI_Scatter and MPI_gather corresponding to gatherList() and scatterList() in Foam.
    http://www.cfd-online.com/Forums/openfoam-solving/58174-reduce-operation-array.html
    */
}

// ************************************************************************* //
