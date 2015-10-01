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
    resLES

Description
    Calculates and writes some mean fields from an LES calculation. The
    following fields are created:
        - resLES: resolutness of an LES
        - TKEMean: mean turbulent kinectic energy
        - TKEMeanProd: production term of the mean turbulent kinectic energy
    Fields UMean, UPrime2Mean and kMean must exist.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    IOobject UPrime2MeanHeader
    (
        "UPrime2Mean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject kMeanHeader
    (
        "kMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject UMeanHeader
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    
  
    if
    (
        UPrime2MeanHeader.headerOk()
     && kMeanHeader.headerOk()
     && UMeanHeader.headerOk()
    )
    {
        Info<< "    Reading average field UPrime2Mean" << endl;
        const volSymmTensorField UPrime2Mean(UPrime2MeanHeader, mesh);

        Info<< "    Reading average field kMean" << endl;
        const volScalarField kMean(kMeanHeader, mesh);

        Info<< "    Reading average field UMean" << endl;
        const volVectorField UMean(UMeanHeader, mesh);


        volScalarField TKEMean
        (
            IOobject
            (
                "TKEMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dimkMean",kMean.dimensions(),0)
        );

        volScalarField resLES
        (
            IOobject
            (
                "resLES",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dimless",dimless,0)
        );

        volScalarField TKEMeanProd
        (
            IOobject
            (
                "TKEMeanProd",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dimTKEMeanProd",kMean.dimensions()/dimTime,0)
        );

        Info<< "    Calculating mean turbulent kinetic energy TKEMean" << endl;
        TKEMean = 0.5 * (
                            UPrime2Mean.component(0)
                          + UPrime2Mean.component(3)
                          + UPrime2Mean.component(5)
                        );

        Info<< "    Calculating Mean LES resolutness resLES" << endl;
        forAll(resLES, i)
        {
            resLES[i] = TKEMean[i] / (TKEMean[i] + kMean[i]);
        }

        Info<< "    Calculating Mean TKE production TKEMeanProd" << endl;
        TKEMeanProd = - UPrime2Mean && fvc::grad(UMean);
        
        resLES.write();
        TKEMean.write();
        TKEMeanProd.write();
    }
    else
    {
        Info << "    No UPrime2Mean and/or kMean and/or UMean." << endl;
    }
    
    Info<< "End" << endl;
}


// ************************************************************************* //
