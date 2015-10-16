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

Application
    TKEP2D

Description
    Compute the Turbulent Kinetic Energy (TKE) in 2D. One direction
    is removed.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    Foam::argList::addOption
    (
        "dir",
        "word",
        "Direction to remove in TKE. -dir x, or y, or z."
    );
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    
    if (args.optionFound("dir")==true)
    {
        word dir(args.optionRead<word>("dir"));
        
        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);
            Info<< "Time = " << runTime.timeName() << endl;
            
            
            IOobject UPrime2MeanHeader
            (
                "UPrime2Mean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            );
            
            volScalarField TKE2D
            (
                IOobject
                (
                    "TKE2D",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "dimTKE2D",
                    pow(dimLength,2)/pow(dimTime,2),
                    0
                ),
                "zeroGradient"
            );
            
            if (UPrime2MeanHeader.headerOk())
            {
                Info<< "    Reading average field UPrime2Mean" << endl;
                volSymmTensorField UPrime2Mean(UPrime2MeanHeader, mesh);         
                
                if (dir=="x")
                {
                    Info << "    Calculating 2D turbulent kinetic energy TKE2D without the " 
                         << dir << "component." << endl;
                    TKE2D = 0.5 * (
//                               UPrime2Mean.component(0)
                              UPrime2Mean.component(3)
                            + UPrime2Mean.component(5)
                            );
                }
                else if (dir=="y")
                {
                    Info << "    Calculating 2D turbulent kinetic energy TKE2D without the " 
                         << dir << "component." << endl;
                    TKE2D = 0.5 * (
                              UPrime2Mean.component(0)
//                             + UPrime2Mean.component(3)
                            + UPrime2Mean.component(5)
                            );
                }
                else if (dir=="z")
                {
                    Info << "    Calculating 2D turbulent kinetic energy TKE2D without the " 
                         << dir << "component." << endl;
                    TKE2D = 0.5 * (
                              UPrime2Mean.component(0)
                            + UPrime2Mean.component(3)
//                             + UPrime2Mean.component(5)
                            );
                }
                else
                {
                    Info << "    dir must be x, y or z" << endl;
                }
            
            TKE2D.write();


            }
            else
            {
                Info << "No UPrime2Mean." << endl;
                return 1;
            }

        } 
    }
    else
    {
        Info << "The 'optional' flag dir must be specified. "
             << "Therefore it is not optional..." << endl;
        return 1;
    }

    

    Info<< nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
