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
    LESdelta

Description
    Calculates and write the LES filter width "LESdelta" (scalar field).

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"
#include "fvCFD.H"
#include "LESModel.H"
#include "singlePhaseTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    Info<< "Reading velocity field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading/calculating face flux field phi\n" << endl;
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        linearInterpolate(U) & mesh.Sf()
    );

    Info<< "Creating LES filter width field LESdelta\n" << endl;
    volScalarField LESdelta
    (
        IOobject
        (
            "LESdelta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("LESdelta", dimLength, SMALL)
    );

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::LESModel> sgsModel
    (
        incompressible::LESModel::New(U, phi, laminarTransport)
    );

    LESdelta.internalField() = sgsModel->delta();

    LESdelta.write();

    Info<< "End" << endl;
}


// ************************************************************************* //

