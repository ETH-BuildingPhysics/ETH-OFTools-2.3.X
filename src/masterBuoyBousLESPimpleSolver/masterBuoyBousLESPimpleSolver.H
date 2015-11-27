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
    Foam::masterBuoyBousLESPimpleSolver

Description
    

SourceFiles
    masterBuoyBousLESPimpleSolver.C

\*---------------------------------------------------------------------------*/

#ifndef masterBuoyBousLESPimpleSolver_H
#define masterBuoyBousLESPimpleSolver_H

#include "runTimeSelectionTables.H"

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "radiationModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"

#include "nestedBlendDsBuoyBousLESPimpleSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mySRFModel Declaration
\*---------------------------------------------------------------------------*/

class masterBuoyBousLESPimpleSolver
{

protected:
    //- masterBuoyBousLESPimpleSolver name
    word name_;

    //- Time object
    Time& runTime_;

    //- mesh
    fvMesh mesh_;

    IOdictionary nestingControlDict_;

    uniformDimensionedVectorField g_;

    //- Temperature field [K]
    volScalarField T_;

    volScalarField p_rgh_;

    //- velocity field [m/s]
    volVectorField U_;

    surfaceScalarField phi_;

    functionObjectList functionObjects_;

    singlePhaseTransportModel laminarTransport_;

    // Thermal expansion coefficient [1/K]
    dimensionedScalar beta_;

    // Reference temperature [K]
    dimensionedScalar TRef_;

    // Laminar Prandtl number
    dimensionedScalar Pr_;

    // Turbulent Prandtl number
    dimensionedScalar Prt_;

    autoPtr<incompressible::turbulenceModel> turbulence_;

    // Kinematic density for buoyancy force
    volScalarField rhok_;

    // kinematic turbulent thermal thermal conductivity m2/s
    volScalarField alphat_;
    
    // Create the wall shear stress field
    volSymmTensorField Rwall_;
    
    // Create the wall temperature flux field
    volVectorField qwall_;

    volScalarField gh_;

    surfaceScalarField ghf_;

    volScalarField p_;

    autoPtr<radiation::radiationModel> radiation_;

    dimensionedScalar rhoCpRef_;

    fv::IOoptionList fvOptions_;

    pimpleControl pimple_;

    label pRefCell_;

    scalar pRefValue_;

    scalar cumulativeContErr_;

    //- list of nested case name (children of this)
    wordList nestedCases_;

    //- list of nested solver (pointers)
    List<nestedBlendDsBuoyBousLESPimpleSolver* > nestedSolvers_;

    bool runNested_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        masterBuoyBousLESPimpleSolver(const masterBuoyBousLESPimpleSolver&);

        //- Disallow default bitwise assignment
        void operator=(const masterBuoyBousLESPimpleSolver&);

public:

    //- Runtime type information
//    TypeName("masterBuoyBousLESPimpleSolver");

    // Declare runtime constructor selection table
//     declareRunTimeSelectionTable
//     (
//         autoPtr,
//         masterBuoyBousLESPimpleSolver,
//         dictionary,
//         (
//             word name, Time& runTime
//         ),
//         (name,runTime)
//     );

    // Constructors

        //- Construct from attributes
        masterBuoyBousLESPimpleSolver
        (
            word name,
            Time& runTime,
            bool runNested = true
        );

        //- Destructor
        virtual ~masterBuoyBousLESPimpleSolver();


    // Member functions

        //- set nestedSolver list
        void setNestedSolvers
        (
            List< nestedBlendDsBuoyBousLESPimpleSolver* > nestedSolvers
        );


        //- run the function objects stored in nestingControlDict
        void runFunctionObjects();

        //- solve the pimple loop
        void solve();


    // Access functions

        //- return name
        word& name();

        //- return nestingControl dictionnary
        IOdictionary& nestingControlDict();

        //- return runTime
        Time& runTime();

        //- return mesh
        fvMesh& mesh();

        //- return temperature field
        volScalarField& T();

        //- return U field
        volVectorField& U();

        //- return p_rgh field
        volScalarField& p_rgh();

        //- return p field
        volScalarField& p();

        //- return phi field
        surfaceScalarField& phi();

        //- return function object list
        functionObjectList& functionObjects();

        //- return the list of nested solver
        List< nestedBlendDsBuoyBousLESPimpleSolver* >& nestedSolvers();

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
