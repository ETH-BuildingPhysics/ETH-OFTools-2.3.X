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


Description


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
    IOobject UMeanHeader
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    
  
    if (UMeanHeader.headerOk())
    {
        Info<< "\nReading field UMean" << endl;
        const volVectorField UMean(UMeanHeader, mesh);

        Info<< "Calculating normal velocity gradient snGradU\n" << endl;
        surfaceVectorField calc_snGradUMean = fvc::snGrad(UMean);
        surfaceVectorField::GeometricBoundaryField& boundary_snGradUMean =
                 calc_snGradUMean.boundaryField();

        volVectorField snGradUMean
        (
            IOobject
            (
                "snGradUMean",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedVector("snGradUMean", calc_snGradUMean.dimensions(), vector::zero)
        );

        forAll(snGradUMean.boundaryField(), patchi)
        {
            snGradUMean.boundaryField()[patchi] = boundary_snGradUMean[patchi];
        }

        snGradUMean.write();
    }
    else
    {
        Info << "No field UMean. snGradUMean exit.\n" << endl;
    }
    
    Info<< "End" << endl;
}


// ************************************************************************* //
