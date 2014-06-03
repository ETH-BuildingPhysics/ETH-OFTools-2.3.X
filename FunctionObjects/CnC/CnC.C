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

\*---------------------------------------------------------------------------*/

#include "CnC.H"
#include "Time.H"
#include "dynamicCode.H"
#include "fvMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(CnC, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CnC::CnC
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool
)
:
    name_(name),
	obr_(obr),
    executeCalls_(),
    endCalls_(),
    writeCalls_(),
	restartCalls_(),
	restartRequested_(false),
	maxClockTime_(-1),
	deltaT_(-1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CnC::~CnC()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CnC::read(const dictionary& dict)
{
	dict.readIfPresent("executeCalls", executeCalls_);
	dict.readIfPresent("endCalls", endCalls_);
	dict.readIfPresent("writeCalls", writeCalls_);
	dict.readIfPresent("restartCalls", restartCalls_);

	
	maxClockTime_ = readScalar(dict.lookup("maxClockTime"));

    if (executeCalls_.empty() && endCalls_.empty() && writeCalls_.empty() && restartCalls_.empty())
    {
        WarningIn("Foam::system::read(const dictionary&)")
            << "no executeCalls, endCalls, writeCalls or restartCalls defined."
            << endl;
    }
    else if (!dynamicCode::allowSystemOperations)
    {
        FatalErrorIn("CnC::read(const dictionary&)")
            << "Executing user-supplied system calls is not enabled by "
            << "default because of " << nl
            << "security issues.  If you trust the case you can enable this "
            << "facility by " << nl
            << "adding to the InfoSwitches setting in the system controlDict:"
            << nl << nl
            << "    allowSystemOperations 1" << nl << nl
            << "The system controlDict is either" << nl << nl
            << "    ~/.OpenFOAM/$WM_PROJECT_VERSION/controlDict" << nl << nl
            << "or" << nl << nl
            << "    $WM_PROJECT_DIR/etc/controlDict" << nl << nl
            << exit(FatalError);
    }
}


void Foam::CnC::execute()
{
	//Info << "CnC: execute()" << endl;

	Time& runTime = const_cast<Time&>(obr_.time());

	scalar elapsedTime = runTime.elapsedClockTime();
	//scalar deltaT = runTime.deltaT().value();
	
	if(elapsedTime >= maxClockTime_)
	{
		Info << "CnC: Clock time is larger than " << maxClockTime_ << " s" << endl;
		Info << "CnC: Sopping and writing data" << endl;
		restartRequested_ = true;
		//runTime.stopAt(Time::saWriteNow);
		runTime.writeAndEnd();
	}

	//Info << "current time:" << runTime.value() << endl;
	//Info << "end Time:" << runTime.endTime() << endl;
	//Info << "dT:" << deltaT << endl;


	/*
	if ((runTime.value()+1.5*deltaT)>runTime.endTime().value())
	{
		Info << "ENDING" << endl;
		scalar ed = runTime.endTime().value() +runTime.deltaT().value();
		//scalar t = 20.0;
		Info << "ed " << ed << endl;
		runTime.setEndTime(ed);
		runTime.stopAt(Time::saWriteNow);
	}*/

	forAll(executeCalls_, callI)
	{
		Foam::system(executeCalls_[callI]);
	}
}


void Foam::CnC::end()
{	
	Info << "CnC: end()" << endl;
	Time& runTime = const_cast<Time&>(obr_.time());
	
	Info << "CnC: writing final timestep" << endl;
	runTime.writeNow();

	forAll(endCalls_, callI)
	{
		Foam::system(endCalls_[callI]);
	}
	
	
	//runTime.writeAndEnd();
	//runTime.stopAt(Time::saWriteNow);
	//runTime.writeNow();
}


void Foam::CnC::timeSet()
{
    // Do nothing
}


void Foam::CnC::write()
{
	Info << "CnC: write()" << endl;
	if (restartRequested_)
	{
		if (Pstream::master())
		{
			Info << "CnC: Calling restart script" << endl;
			forAll(restartCalls_, callI)
			{
				Foam::system(restartCalls_[callI]);
			}
		}
		
	}else
	{
		if (Pstream::master())
		{
			Info << "CnC: Calling write script" << endl;
			forAll(writeCalls_, callI)
			{
				Foam::system(writeCalls_[callI]);
			}
		}
	}
}


// ************************************************************************* //
