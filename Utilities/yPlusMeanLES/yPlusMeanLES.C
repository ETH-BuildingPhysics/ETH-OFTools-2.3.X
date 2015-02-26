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
    yPlusMeanLES

Description
    Calculates and reports yPlusMean for all wall patches, for the specified
    times when using LES turbulence models. UMean must exist.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "LESModel.H"
#include "nearWallDist.H"
#include "wallDist.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        fvMesh::readUpdateState state = mesh.readUpdate();

        // Wall distance
        if (timeI == 0 || state != fvMesh::UNCHANGED)
        {
            Info<< "Calculating wall distance\n" << endl;
            wallDist y(mesh, true);
            Info<< "Writing wall distance to field "
                << y.name() << nl << endl;
            y.write();
        }


        volScalarField yPlusMean
        (
            IOobject
            (
                "yPlusMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlusMean", dimless, 0.0)
        );


        Info<< "Reading field UMean\n" << endl;
        volVectorField UMean
        (
            IOobject
            (
                "UMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        Info<< "Reading field nuSgsMean\n" << endl;
        volScalarField nuSgsMean
        (
            IOobject
            (
                "nuSgsMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        #include "createPhiMean.H"

        singlePhaseTransportModel laminarTransport(UMean, phiMean);

        autoPtr<incompressible::LESModel> sgsModel
        (
            incompressible::LESModel::New(UMean, phiMean, laminarTransport)
        );
        const volScalarField nuLam(sgsModel->nu());
        volScalarField nuEff(nuLam+nuSgsMean);

         volScalarField::GeometricBoundaryField d = nearWallDist(mesh).y();
        const fvPatchList& patches = mesh.boundary();


        forAll(patches, patchi)
        {
            const fvPatch& currPatch = patches[patchi];

            if (isA<wallFvPatch>(currPatch))
            {
                yPlusMean.boundaryField()[patchi] =
                    d[patchi]
                   *sqrt
                    (
                        nuEff.boundaryField()[patchi]
                       *mag(UMean.boundaryField()[patchi].snGrad())
                    )
                   /nuLam.boundaryField()[patchi];
                const scalarField& Yp = yPlusMean.boundaryField()[patchi];

                Info<< "Patch " << patchi
                    << " named " << currPatch.name()
                    << " y+ : min: " << gMin(Yp) << " max: " << gMax(Yp)
                    << " average: " << gAverage(Yp) << nl << endl;
            }
        }

        Info<< "Writing yPlusMean to field "
            << yPlusMean.name() << nl << endl;

        yPlusMean.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
