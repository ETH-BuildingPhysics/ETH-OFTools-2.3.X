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
    nestedImpBlendDsBuoyBousLESPimpleFoam

Description
    One-way, downscale, nested solver. A nested (or inner) solver is
    embedded in a master (or outer) solver. The outter solver has a coarser
    grid and timestep whereas the inner domain runs on a finer grid and
    timestep. The inner grid is nested to the outer one with nested boundaries
    and a blending zone. Such procedure is particularly useful to simulate the
    flow around man made structures (small scale) fully immersed in an 
    Atmospheric Boundary Layer (ABL). See tutorial for more details on usage.
    
    The buoyant effects are modeled with the Boussinescq approximation for
    buoyancy. The fields Qwall and rwall are used for LES wall-model. Despite
    beibg designed primary for LES-to-LES nesting, this solver works
    perfectly for URANS-to-URANS nesting. The blending parameters should be 
    adapted (the relaxation time can be higher I guess...)
    
Publication
    Solver based on the following publication:
        Vonlanthen M. et al. Submitted (2016). "Assessment of a One-way Nesting
        Procedure for Obstacle Resolved Large Eddy Simulation of the ABL". 
        Computers and Fluids.

Author
    Marcel Vonlanthen (vonlanthen[dot]marcel[at]gmail[dot]com)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"

#include "masterBuoyBousLESPimpleSolver.H"
#include "nestedBlendDsBuoyBousLESPimpleSolver.H"
#include "NestedTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Create global variable for the run
    #include "setRootCase.H"
    #include "createTime.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    IOdictionary nestingControlDict
    (
         IOobject
         (
              "nestingControlDict",
              runTime.system(),
              runTime,
              IOobject::MUST_READ,
              IOobject::NO_WRITE,
              true
         )
    );

//    // load the master case
    word masterName(Foam::fvMesh::defaultRegion);
    Info << nl << "Load case: " << masterName << nl << endl;
    masterBuoyBousLESPimpleSolver masterSolver(masterName,runTime);


    // load the nested cases
    wordList allNestedCases(nestingControlDict.lookup("allNestedCases"));
    List< nestedBlendDsBuoyBousLESPimpleSolver* > allNestedSolvers(allNestedCases.size());
    forAll (allNestedCases,i)
    {
        word namei = allNestedCases[i];
        Info << nl << "Load case: " << namei << nl << endl;
        allNestedSolvers[i] = new nestedBlendDsBuoyBousLESPimpleSolver(namei,runTime,true);

        // normally, the methods "creatInterpolator()" and "setwDs()" of the
        // nestedSolver should always be executed. Nevertheless meshToMesh
        // has a bug.
        IOdictionary& ntdCtlDicti(allNestedSolvers[i]->nestingControlDict());
        bool scaleMesh(ntdCtlDicti.lookupOrDefault("scaleMesh",false));
        bool doBlending(ntdCtlDicti.lookupOrDefault("doBlending",true));
        if (scaleMesh==false && doBlending==true)
        {
            allNestedSolvers[i]->creatInterpolator();
            allNestedSolvers[i]->setwDs();
        }
        if (scaleMesh==true && doBlending==true)
        {
            scalar scaleFactor(ntdCtlDicti.lookupOrDefault("scaleFactor",0.001));

            // get the raw and the scaled master points
            const pointField& rawPtMst = masterSolver.mesh().points();
            pointField orgPtMst = masterSolver.mesh().points();

            // get the raw and the scaled nested points
            const pointField& rawPtNtd =  allNestedSolvers[i]->mesh().points();
            pointField orgPtNtd =  allNestedSolvers[i]->mesh().points();

            // scale down the points
            Info << "scale meshes down!" << endl;
            masterSolver.mesh().movePoints(scaleFactor*rawPtMst);
            allNestedSolvers[i]->mesh().movePoints(scaleFactor*rawPtNtd);

            // do the adressing
            allNestedSolvers[i]->creatInterpolator();

            // scale up the points
            Info << "scale meshes up!" << endl;
            masterSolver.mesh().movePoints(orgPtMst);
            allNestedSolvers[i]->mesh().movePoints(orgPtNtd);

            masterSolver.mesh().moving(false);
            allNestedSolvers[i]->mesh().moving(false);

            // set wDs
            allNestedSolvers[i]->setwDs();
        }
    };

    // set the list of nested solver in each solver
    masterSolver.setNestedSolvers(allNestedSolvers);
    forAll (allNestedSolvers,i)
    {
        allNestedSolvers[i]->setNestedSolvers(allNestedSolvers);
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        runTime++;

        // solve pimple loop of the master case
        masterSolver.solve();

        runTime.write();
    };


    Info<< "End\n" << endl;

    return 0;

}
// ************************************************************************* //
