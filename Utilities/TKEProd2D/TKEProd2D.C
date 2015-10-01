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
    TKEProd2D

Description
    Compute the Turbulent Kinetic Energy (TKE) production in 2D. One direction
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
        "Direction to remove in the TKE production term. -dir x, or y, or z."
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
            
            IOobject UMeanHeader
            (
                "UMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ  
            );
            
            IOobject UPrime2MeanHeader
            (
                "UPrime2Mean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            );
            
            if
            (
                UMeanHeader.headerOk()
            && UPrime2MeanHeader.headerOk()
            )
            {
                Info<< "    Reading average field UMean" << endl;
                const volVectorField UMean(UMeanHeader, mesh);
                
                Info<< "    Reading average field UPrime2Mean" << endl;
                volSymmTensorField UPrime2Mean(UPrime2MeanHeader, mesh);
                
                volTensorField gradUMean(fvc::grad(UMean));          
                
                if (dir=="x")
                {
                    gradUMean.replace(0,pTraits<scalar>::zero);
                    gradUMean.replace(3,pTraits<scalar>::zero);
                    gradUMean.replace(6,pTraits<scalar>::zero);
                    gradUMean.replace(1,pTraits<scalar>::zero);
                    gradUMean.replace(2,pTraits<scalar>::zero);

                    UPrime2Mean.replace(0,pTraits<scalar>::zero);
                    UPrime2Mean.replace(1,pTraits<scalar>::zero);
                    UPrime2Mean.replace(2,pTraits<scalar>::zero);
                }
                else if (dir=="y")
                {
                    gradUMean.replace(1,pTraits<scalar>::zero);
                    gradUMean.replace(4,pTraits<scalar>::zero);
                    gradUMean.replace(7,pTraits<scalar>::zero);
                    gradUMean.replace(3,pTraits<scalar>::zero);
                    gradUMean.replace(5,pTraits<scalar>::zero);

                    UPrime2Mean.replace(1,pTraits<scalar>::zero);
                    UPrime2Mean.replace(3,pTraits<scalar>::zero);
                    UPrime2Mean.replace(4,pTraits<scalar>::zero);
                }
                else if (dir=="z")
                {
                    gradUMean.replace(2,pTraits<scalar>::zero);
                    gradUMean.replace(5,pTraits<scalar>::zero);
                    gradUMean.replace(8,pTraits<scalar>::zero);
                    gradUMean.replace(6,pTraits<scalar>::zero);
                    gradUMean.replace(7,pTraits<scalar>::zero);

                    UPrime2Mean.replace(3,pTraits<scalar>::zero);
                    UPrime2Mean.replace(4,pTraits<scalar>::zero);
                    UPrime2Mean.replace(5,pTraits<scalar>::zero);
                }
                else
                {
                    Info << "    dir must be x, y or z" << endl;
                }

                volScalarField TKEProd2D
                (
                    IOobject
                    (
                        "TKEProd2D",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar
                    (
                        "dimTKEProd2D",
                        pow(dimLength,2)/pow(dimTime,3),
                        0
                    ),
                    "zeroGradient"
                );
                TKEProd2D = - UPrime2Mean && fvc::grad(UMean);
                Info << "    Write field TKEProd2D with '" << dir
                     <<      "' the removed direction." << endl;
                TKEProd2D.write();
            }
            else
            {
                Info << "No UMean and/or UPrime2Mean." << endl;
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
