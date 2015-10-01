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

\*---------------------------------------------------------------------------*/

#include "NestedTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
//        defineTypeNameAndDebug(NestedTime, 0);
//        defineRunTimeSelectionTable(NestedTime, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

NestedTime::NestedTime
(
    const word& name,
    const argList& args
//    const word& systemName = "system",
//    const word& constantName = "constant"
)
:
    Time
    (
        name,
        args
//        systemName,
//        constantName
    ),
    LIFOprevTimeState_(),
    LIFOdeltaT_(),
    LIFOdeltaT0_(),
    LIFOdeltaTSave_()
{
    LIFOprevTimeState_.begin();
    LIFOdeltaT_.begin();
    LIFOdeltaT0_.begin();
    LIFOdeltaTSave_.begin();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

NestedTime::~NestedTime()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::TimeState Foam::NestedTime::subCycle(const label nSubCycles)
{
    subCycling_ = true;
    prevTimeState_.clear();
    prevTimeState_.set(new TimeState(*this));
    LIFOprevTimeState_.push( prevTimeState_ );

    setTime(*this - deltaT(), (timeIndex() - 1)*nSubCycles);
    deltaT_ /= nSubCycles;
    deltaT0_ /= nSubCycles;
    deltaTSave_ = deltaT0_;

    LIFOdeltaT_.push(deltaT_);
    LIFOdeltaT0_.push(deltaT0_);
    LIFOdeltaTSave_.push(deltaTSave_);

    return prevTimeState();
}


void Foam::NestedTime::endSubCycle()
{
    if (subCycling_)
    {
        TimeState::operator=(LIFOprevTimeState_.pop());
        if (LIFOprevTimeState_.size()==0)
        {
            subCycling_ = false;
            prevTimeState_.clear();
        }

        LIFOdeltaT_.pop();
        LIFOdeltaT0_.pop();
        LIFOdeltaTSave_.pop();
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
