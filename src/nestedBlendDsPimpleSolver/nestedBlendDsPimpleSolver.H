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
    Foam::nestedBlendDsPimpleSolver

Description
    A class, which handles nested solver of a nested one-way downscale
    LES-within-LES solver. See class masterPimpleSolver for the master part 
    of the nested procedure. This class is mainly used to construct a nested
    solver. The turbulence model can be laminar, RANS or LES.
    
SourceFiles
    nestedBlendDsPimpleSolver.C
    
Author
    Marcel Vonlanthen (vonlanthen[dot]marcel[at]gmail[dot]com)

\*---------------------------------------------------------------------------*/

#ifndef nestedBlendDsPimpleSolver_H
#define nestedBlendDsPimpleSolver_H

#include "runTimeSelectionTables.H"

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"

#include "subCycleTime.H"
#include "patchWave.H"

#include "nestingDsInterpolator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class nestedBlendDsPimpleSolver Declaration
\*---------------------------------------------------------------------------*/

class nestedBlendDsPimpleSolver
{
    // main member variables
    
        //- name of the solver. Passed to the mesh in the database
        word name_;

        //- reference to the time object
        Time& runTime_;

        //- mesh
        fvMesh mesh_;

        // the nestingControlDict parser
        IOdictionary nestingControlDict_;

        //- pressure field
        volScalarField p_;

        //- velocity field
        volVectorField U_;

        //- momentum flux
        surfaceScalarField phi_;

        //- list of functionObjects defined in nestingControlDict
        functionObjectList functionObjects_;

        //- the transport model
        singlePhaseTransportModel laminarTransport_;

        //- auto-pointer of the turbulent model
        autoPtr<incompressible::turbulenceModel> turbulence_;

        //- fvOptions
        fv::IOoptionList fvOptions_;

        //- pimple control
        pimpleControl pimple_;

        //- reference cell and pressure for the pimple solver
        label pRefCell_;

        scalar pRefValue_;

        //- cumulative error
        scalar cumulativeContErr_;

    // data from the nested cases (children of this)

        //- list of nested case name.
        wordList nestedCases_;

        //- list of nested solver (pointers)
        List< nestedBlendDsPimpleSolver* > nestedSolvers_;

        bool runNested_;

    // data from the master

        //- master name
        const word nameMst_;

        //- master mesh
        const fvMesh& meshMst_;

        //- downscale interpolator
        nestingDsInterpolator interpolator_;

    // extra field for downscale nesting
        
        //- velocity field of the master interpolated on
        //  this mesh
        volVectorField UMstOnNtd_;

        volVectorField U0MstOnNtd_;

        volVectorField UDsBlended_;

        volScalarField pMstOnNtd_;

        volScalarField p0MstOnNtd_;

        volScalarField pDsBlended_;

        //- blending factor
        volScalarField wDs_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        nestedBlendDsPimpleSolver(const nestedBlendDsPimpleSolver&);

        //- Disallow default bitwise assignment
        void operator=(const nestedBlendDsPimpleSolver&);

public:

    //- Runtime type information
//    TypeName("nestedBlendDsPimpleSolver");

    // Declare runtime constructor selection table
//     declareRunTimeSelectionTable
//     (
//         autoPtr,
//         nestedBlendDsPimpleSolver,
//         dictionary,
//         (
//             const fvMesh& mesh
//         ),
//         (mesh)
//     );

    // Constructors

        //- Construct from attributes
        nestedBlendDsPimpleSolver
        (
            word name,
            Time& runTime,
            bool runNested = true
        );

        //- Destructor
        virtual ~nestedBlendDsPimpleSolver();


    // Member functions

        //- update wDs field according
        void setwDs();

        //- interlpolate field
        void interpolateMstFields();

        //- set nestedSolver list
        void setNestedSolvers
        (
            List< nestedBlendDsPimpleSolver* > nestedSolvers
        );

        //- run the function objects stored in nestingControlDict
        void runFunctionObjects();

        //- solve the pimple loop
        void solve();


    // Access functions

        //- return name
        word& name();

        //- return nestingControl
        IOdictionary& nestingControlDict();

        //- return runTime
        Time& runTime();

        //- return mesh
        fvMesh& mesh();

        //- return U field
        volVectorField& U();

        //- return p field
        volScalarField& p();

        //- return phi field
        surfaceScalarField& phi();

        //- return function object list
        functionObjectList& functionObjects();

        //- return
        List< nestedBlendDsPimpleSolver* >& nestedSolvers();

        //- return if the nested solvers are run
        bool runNested();


    // setter functions

        //- set runNested
        void runNested(bool b);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
