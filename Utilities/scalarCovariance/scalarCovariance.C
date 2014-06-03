/*---------------------------------------------------------------------------* \
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
    Calculates and writes <u's'> = <us> - <u><s>

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
	bool writeResults = !args.optionFound("noWrite");

	const word EXT_MEAN = "Mean";
	const word EXT_PRIME = "Prime";
	const word EXT_PRIME2MEAN = "Prime2Mean";

	#include "createFields.H"

	forAll(names, i)
        {
		Info<< "    Calculating "+names[i] << endl;
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

		if (writeResults)
		{
		    UsPrime.write();
		}
	}

	Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
