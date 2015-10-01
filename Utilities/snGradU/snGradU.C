/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    snGradP

Description
    Calculates and writes the normal pressure gradient snGradP for each patches.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"
#include "fvCFD.H"
// #include "dimensionedVector.H"
// #include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
// #include "RASModel.H"
// #include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    IOobject UHeader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    
  
    if (UHeader.headerOk())
    {
        Info<< "\nReading field U" << endl;
        const volVectorField U(UHeader, mesh);

        Info<< "Calculating normal velocity gradient snGradU\n" << endl;
        surfaceVectorField calc_snGradU = fvc::snGrad(U);
        surfaceVectorField::GeometricBoundaryField& boundary_snGradU =
                 calc_snGradU.boundaryField();

        volVectorField snGradU
        (
            IOobject
            (
                "snGradU",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedVector("snGradU", calc_snGradU.dimensions(), vector::zero)
        );

        forAll(snGradU.boundaryField(), patchi)
        {
            snGradU.boundaryField()[patchi] = boundary_snGradU[patchi];
        }

        snGradU.write();
    }
    else
    {
        Info << "No field U. snGradU exit.\n" << endl;
    }
    
    Info<< "End" << endl;
}


// ************************************************************************* //
