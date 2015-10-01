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
    scalarCovariance

Description
    For each scalar field Si, Calculates and writes <u'si'> = <USi> - <U><Si>.
    Si can be any volScalarField: passive tracer, T, k, nuSgs...

    The list of scalar field Si can be specified in transportProperty as a list:
    scalars (S1 S2 ... Si ... Sn), Or as an argument:
    scalarConvariance -scalars "(S1 S2 ... Si ... Sn)".

    For any scalar field Si, the following fields must exist: UMean, SiMean,
    USiMean. All those fields can be generated in run-time with the function
    objects "scalarVelocityProduct" (for USi) and "fieldAverage"
    (for the means).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    Foam::argList::addOption
    (
        "scalars",
        "list",
        "specify a list of scalar. \"(S1..Si..Sn)\" or Si for a single scalar."
    );
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    const word EXT_MEAN = "Mean";
    const word EXT_PRIME = "Prime";
    const word EXT_PRIME2MEAN = "Prime2Mean";

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        #include "createFields.H"

        forAll(names, i)
        {
            Info << "Calculating field "
                 << "U"+EXT_PRIME+names[i]+EXT_PRIME+EXT_MEAN << endl;

            volVectorField UsPrime
            (
                IOobject
                (
                    "U"+EXT_PRIME+names[i]+EXT_PRIME+EXT_MEAN,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                usMean[i] - (UMean * sMean[i])
            );

            UsPrime.write();
        }
    }

    Info<< nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
