/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    patchDistance

Description
    Compute the distance from a patch or a list of patches and save it
    in the field 'yDist_nameExtension'. The list of patches and the
    nameExtension are defined in the dictonnary 'patchDistanceDict' in
    'system'.
    
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "patchWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    
    IOdictionary patchDistanceDict
    (
        IOobject
        (
            "patchDistanceDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    const wordReList srcPatches(patchDistanceDict.lookup("sourcePatches"));
    word fieldNameExt(patchDistanceDict.lookup("fieldNameExtension"));
    word fieldName = "yDist_"+fieldNameExt;
    
    Info << "yDist computed from patches:" << endl;
    forAll(srcPatches, patchI)
    {
        Info << "   " << srcPatches[patchI] << endl;
    };
    Info << "and saved with the field name '" << fieldName << "'\n" << endl;
    
    #include "createNamedMesh.H"
    
    labelHashSet srcPatchesSet(mesh.boundaryMesh().patchSet(srcPatches));
    patchWave yDistWave(mesh,srcPatchesSet,true);
    volScalarField yDist
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("dimLength",dimLength,0),
        "zeroGradient"                         // zeroGradient at the boundary
    );
    yDist.internalField() = yDistWave.distance();
    Info << "write " << fieldName << endl;
    yDist.write();

    Info<< nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
