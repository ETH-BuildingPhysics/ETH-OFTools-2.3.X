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
    Foam::NestedTime

Description
    Derived class of runTime. It can handle multiple sub-cycles embbeded
    into each other (a sub-sub-cycle). This class is used in the nested solvers.
    It can be used exctly like the Foam::Time class.
    
SourceFiles
    NestedTime.C
    
Author
    Marcel Vonlanthen (vonlanthen[dot]marcel[at]gmail[dot]com)

\*---------------------------------------------------------------------------*/

#ifndef NestedTime_H
#define NestedTime_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "runTimeSelectionTables.H"
#include "meshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mySRFModel Declaration
\*---------------------------------------------------------------------------*/

class NestedTime
:
    public Time
{
    // Private data


protected:
    // Protected data

        //- LIFO stack of previous time state. Used for multiple subCycle
        LIFOStack< TimeState > LIFOprevTimeState_;

        //- LIFO stack of deltaT
        LIFOStack<scalar> LIFOdeltaT_;

        //- LIFO stack of the previous deltaT
        LIFOStack<scalar> LIFOdeltaT0_;
        
        //- 
        LIFOStack<scalar> LIFOdeltaTSave_;




    //- IO dictonnary
//    IOdictionary nestingControlDict_;

    //- master mesh
//    const fvMesh& meshMst_;



private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        NestedTime(const NestedTime&);

        //- Disallow default bitwise assignment
        void operator=(const NestedTime&);


public:

    //- Runtime type information
//    TypeName("NestedTime");

    // Declare runtime constructor selection table
//     declareRunTimeSelectionTable
//     (
//         autoPtr,
//         NestedTime,
//         dictionary,
//         (
//             const fvMesh& meshMst,
//             const fvMesh& meshNtd,
//             const IOdictionary& nestingControlDict
//         ),
//         (meshMst, meshNtd, nestingControlDict)
//     );

    // Constructors

        // Construct from attributes
        NestedTime
        (
            const word& name,
            const argList& args
//            const word& systemName = "system",
//            const word& constantName = "constant"
        );

        //- Destructor
        virtual ~NestedTime();


    // Member functions

        //- Set time to sub-cycle for the given number of steps
        virtual TimeState subCycle(const label nSubCycles);

        //- Reset time after sub-cycling back to previous TimeState
        virtual void endSubCycle();

    // Access functions

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
