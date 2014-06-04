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

Class
    Foam::timeAbort

Group
    grpJobControlFunctionObjects
Author
    Marcel Vonlanthen
    vonlanthen.marcel@outlook.com
    17.03.2014

Description
    Writes and stoppes the solver when the time type (timeType)
    is bigger then the prescribed time value (timeValue or timeVector).
    Useful if OpenFOAM runs with time limited queues, it allows to stop and
    write just before been killed by the queuing system.

    Usage (at the end of the controlDict):
    functions
    {
        timeAbort
        {
            type timeAbort;
            functionObjectLibs ( "libtimeAbort.so" );
            timeType clockTime; // clockTime (default) or executionTime
            // Select one of the two possibilities: timeValue or timeVector.
            // If both are selected, timeValue ONLY is used.
            timeVector (2 35 30);      // time from solver start (h min sec)
            timeValue 240;     //time from solver start (sec)
        }
    }

    Currently the following time types are supported:
    - clockTime
    - executionTime

SourceFiles
    timeAbort.C
    IOtimeAbort.H

\*---------------------------------------------------------------------------*/

#ifndef timeAbort_H
#define timeAbort_H

#include "NamedEnum.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                      Class timeAbort Declaration
\*---------------------------------------------------------------------------*/

class timeAbort
{
public:

    // Public data

        //- emuneration defining the type of time
        enum timeType
        {
            clockTime,
            executionTime
        };

private:

    // Private data

        //- Name of the abort file unless otherwise specified. No clue about
        // this shit!! It's present in all func obj with always a different
        // comment...
        word name_;

        const objectRegistry& obr_;

        //- time type names. two choices: clockTime and executionTime
        static const NamedEnum<timeType, 2> timeTypeNames_;

        //- The type of action
        timeType timeType_;

        //- time Vector which defines when the job must stop (hour min sec)
        vector timeVector_;

        //- time value in second. Caluclated from timeVector_
        scalar timeValue_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        timeAbort(const timeAbort&);

        //- Disallow default bitwise assignment
        void operator=(const timeAbort&);


public:

    //- Runtime type information
    TypeName("timeAbort");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        timeAbort
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFilesUnused = false
        );


    //- Destructor
    virtual ~timeAbort();


    // Member Functions

        //- Return name of the abort file
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the dictionary settings
        virtual void read(const dictionary&);

        //- Execute, check existence of abort file and take action
        virtual void execute();

        //- Execute at the final time-loop, used for cleanup
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Execute, check existence of abort file and take action
        virtual void write();

        //- Update for changes of mesh - does nothing
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh - does nothing
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
