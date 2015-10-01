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

\*---------------------------------------------------------------------------*/

#include "nestingDsInterpolator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
        defineTypeNameAndDebug(nestingDsInterpolator, 0);
        defineRunTimeSelectionTable(nestingDsInterpolator, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nestingDsInterpolator::nestingDsInterpolator
(
    const fvMesh& meshMst,
    const fvMesh& meshNtd,
    const IOdictionary& nestingControlDict
)
:
    IOdictionary
    (
        IOobject
        (
            "nestingDsInterpolator",
            meshNtd.time().system(),
            meshNtd,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    nestingControlDict_(nestingControlDict),
    meshMst_(meshMst),
    meshNtd_(meshNtd),
    mapMethod_(meshToMesh::imCellVolumeWeight),
    patchMap_(nestingControlDict_.lookup("patchMap")),
    cuttingPatches_(nestingControlDict_.lookup("cuttingPatches"))
{
//    lookup("patchMap") >> patchMap_;
//    lookup("cuttingPatches") >>  cuttingPatches_;

//    DSinterpolator_.set
//    (
//        new meshToMesh
//        (
//            meshMst_,
//            meshNtd_,
//            mapMethod_,
//            patchMap_,
//            addProcessorPatches()

//        )
//    );

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

nestingDsInterpolator::~nestingDsInterpolator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void nestingDsInterpolator::creatInterpolator()
{
    DSinterpolator_.set
    (
        new meshToMesh
        (
            meshMst_,
            meshNtd_,
            mapMethod_,
            patchMap_,
            addProcessorPatches()

        )
    );
}



wordList nestingDsInterpolator::addProcessorPatches()
{
    // Add the processor patches to the cutting list
    HashSet<word> cuttingPatchTable;
    forAll(cuttingPatches_, i)
    {
        cuttingPatchTable.insert(cuttingPatches_[i]);
    }

    const polyBoundaryMesh& pbm = meshNtd_.boundaryMesh();

    forAll(pbm, patchI)
    {
        if (isA<processorPolyPatch>(pbm[patchI]))
        {
            const word& patchName = pbm[patchI].name();
            cuttingPatchTable.insert(patchName);
        }
    }

    return cuttingPatchTable.toc();
}

const meshToMesh& nestingDsInterpolator::DSinterpolatorObj() const
{
    return DSinterpolator_;
}

const fvMesh& nestingDsInterpolator::meshMst() const
{
    return meshMst_;
}

const fvMesh& nestingDsInterpolator::meshNtd() const
{
    return meshNtd_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
