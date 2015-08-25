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

Class
    Foam::nestingDsInterpolator

Description
    Class to handle the interpolation process of nested solvers with
    blending zones.

SourceFiles
    nestingDsInterpolator.C

\*---------------------------------------------------------------------------*/

#ifndef nestingDsInterpolator_H
#define nestingDsInterpolator_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "runTimeSelectionTables.H"
#include "meshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nestingDsInterpolator Declaration
\*---------------------------------------------------------------------------*/

class nestingDsInterpolator
:
    public IOdictionary
{
    // Private data

    //- IO dictonnary
    IOdictionary nestingControlDict_;

    //- master mesh
    const fvMesh& meshMst_;

    //- nested mesh
    const fvMesh& meshNtd_;

    //- map method for the mest to mesh interpolation
    const meshToMesh::interpolationMethod mapMethod_;

    //- patch map: pair of patches (master patch / nested patch)
    //  with the same topology.
    HashTable<word> patchMap_;

    //- list of patches of the nested domain which cut the master domain
    wordList cuttingPatches_;

    // pointer to a meshToMesh object
    autoPtr<meshToMesh> DSinterpolator_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        nestingDsInterpolator(const nestingDsInterpolator&);

        //- Disallow default bitwise assignment
        void operator=(const nestingDsInterpolator&);


public:

    //- Runtime type information
    TypeName("nestingDsInterpolator");

    // Declare runtime constructor selection table
     declareRunTimeSelectionTable
     (
         autoPtr,
         nestingDsInterpolator,
         dictionary,
         (
             const fvMesh& meshMst,
             const fvMesh& meshNtd,
             const IOdictionary& nestingControlDict
         ),
         (meshMst, meshNtd, nestingControlDict)
     );

    // Constructors

        // Construct from attributes
        nestingDsInterpolator
        (
            const fvMesh& meshMst,
            const fvMesh& meshNtd,
            const IOdictionary& nestingControlDict
        );

        //- Destructor
        virtual ~nestingDsInterpolator();


    // Member functions

        //- create meshtoMesh object (this is a costly operation).
        //  Normally, this method should be called in the constructor.
        //  But the underlaying meshToMesh class has some (minor) bug, therefore
        //  this method is separated.
        void creatInterpolator();

        //- add processor patches to cuttingPatches_
        wordList addProcessorPatches();


        // Access functions

            //- return interpolator object
            const meshToMesh& DSinterpolatorObj() const;

            //- return master mesh
            const fvMesh& meshMst() const;

            //- return nested mesh
            const fvMesh& meshNtd() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
