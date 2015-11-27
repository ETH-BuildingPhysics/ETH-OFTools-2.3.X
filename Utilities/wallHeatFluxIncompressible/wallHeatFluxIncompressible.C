/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    wallHeatFluxIncompressible

Description
    Calculates and writes the heat flux for all patches as the boundary field
    of a volScalarField and also prints the integrated flux for all wall
    patches.
    Based on wallHeatFlux with changes to allow it on incompressible flows
    Also removed a bug at the typeid checkline
    Eelco van Vliet
    
    
    
    Modification by Marcel Vonlanthen
    I made this utility lighter in term of outputs. gradT is not anymore 
    calculate. I also remove the writing of gradT, which was impossible to plot
    in paraview

    note: The wallheatflux is calculated as:

    wallHeatFlux = alphaEff*rho*Cp*grad(T) 
    
    with alphaEff = nu/Pr + nut/Prt

    Prt, the turbulent prantld number, is assumed to be constant.
\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
//     timeSelector::addOptions();
//     #include "addRegionOption.H"
//     #include "setRootCase.H"
//     #include "createTime.H"
//     instantList timeDirs = timeSelector::select0(runTime, args);
//     #include "createMesh.H"

//     forAll(timeDirs, timeI)
//     {
//         runTime.setTime(timeDirs[timeI], timeI);
//         Info<< "Time = " << runTime.timeName() << endl;
//         mesh.readUpdate();

    #include "createFields.H"
    #include "readTransportProperties.H"

    // update the turbulence fields
    turbulence->read();

    if (!(IOobject("alphat", runTime.timeName(), mesh).headerOk()))
    {
        Info<< "\nCalculating turbulent heat conductivity " << endl;
        alphat = turbulence->nut()/Prt;
        alphat.correctBoundaryConditions();
    }
    else
    {
        Info<< "\nRead turbulent heat conductivity alphat" << endl;
    }

    if (!(IOobject("alphaEff", runTime.timeName(), mesh).headerOk()))
    {
        Info<< "\nCalculating effective heat conductivity " << endl;
        alphaEff=turbulence->nu()/Pr+alphat;
    }
    else
    {
        Info<< "\nRead effective heat conductivity alphaEff" << endl;
    }

    // calculate heatFlux at the surface
    surfaceScalarField heatFlux =fvc::interpolate(alphaEff*Cp*rho)*fvc::snGrad(T);
        
    const surfaceScalarField::GeometricBoundaryField& patchHeatFlux =
                heatFlux.boundaryField();

    // display sum of heat flux
    Info<< "\nWall heat fluxes " << endl;
    forAll(patchHeatFlux, patchi)
    {
        if (typeid(mesh.boundary()[patchi]) == typeid(wallFvPatch))
        {
            Info<< mesh.boundary()[patchi].name()
                << ": Total "
                << sum
                    (
                        mesh.magSf().boundaryField()[patchi]
                        *patchHeatFlux[patchi]
                    )
                << " [W] over "
                << sum
                    (
                        mesh.magSf().boundaryField()[patchi]
                    )
                << " [m2] ("
                << sum
                    (
                        mesh.magSf().boundaryField()[patchi]
                        *patchHeatFlux[patchi]
                    )/
                    sum 
                    (
                        mesh.magSf().boundaryField()[patchi]
                    )
                << " [W/m2])"
                << endl;
        }
    }
    Info<< endl;

    // create wallHeatFlux
    volScalarField wallHeatFlux
    (
        IOobject
        (
            "wallHeatFlux",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("wallHeatFlux", heatFlux.dimensions(), 0.0)
    );

    forAll(wallHeatFlux.boundaryField(), patchi)
    {
        wallHeatFlux.boundaryField()[patchi] = patchHeatFlux[patchi];
    }

    wallHeatFlux.write();
//     }

    Info<< "End" << endl;
}

// ************************************************************************* //
