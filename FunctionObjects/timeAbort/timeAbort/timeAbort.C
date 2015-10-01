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

#include "timeAbort.H"
#include "dictionary.H"
#include "error.H"
#include "Time.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(timeAbort, 0);
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::timeAbort::timeType,
        2
    >::names[] =
    {
        "clockTime",
        "executionTime",
    };
}

const Foam::NamedEnum<Foam::timeAbort::timeType, 2>
    Foam::timeAbort::timeTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeAbort::timeAbort
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    timeType_(clockTime)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeAbort::~timeAbort()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeAbort::read(const dictionary& dict)
{
    if (dict.found("timeVector") && dict.found("timeValue"))
    {
        timeValue_ = readScalar(dict.lookup("timeValue"));
        timeVector_ = vector::zero;
    }
    else if (dict.found("timeVector")==false && dict.found("timeValue"))
    {
        timeValue_ = readScalar(dict.lookup("timeValue"));
        timeVector_ = vector::zero;
    }
    else
    {
        timeVector_ = dict.lookup("timeVector");
        timeValue_ = 3600*timeVector_[0]
                   + 60*timeVector_[1]
                   + timeVector_[2];
    }

    if (dict.found("timeType"))
    {
        timeType_ = timeTypeNames_.read(dict.lookup("timeType"));
    }
    else
    {
        timeType_ = clockTime;
    }
}


void Foam::timeAbort::execute()
{
    switch (timeType_)
    {
        case clockTime :
        {
            if (obr_.time().elapsedClockTime() > timeValue_)
            {
                if (obr_.time().stopAt(Time::saWriteNow))
                {
                    Info<< "USER REQUESTED ABORT (timeIndex="
                        << obr_.time().timeIndex()
                        << "): stop+write data"
                        << endl;
                }
                break;
            }
        }

        case executionTime :
        {
            if (obr_.time().elapsedCpuTime() > timeValue_)
            {
                if (obr_.time().stopAt(Time::saWriteNow))
                {
                    Info<< "USER REQUESTED ABORT (timeIndex="
                        << obr_.time().timeIndex()
                        << "): stop+write data"
                        << endl;
                }
                break;
            }
        }
    }
}


void Foam::timeAbort::end()
{
    // Do nothing - only valid on execute
}


void Foam::timeAbort::timeSet()
{
    // Do nothing - only valid on execute
}


void Foam::timeAbort::write()
{
    // Do nothing - only valid on execute
}


// ************************************************************************* //
