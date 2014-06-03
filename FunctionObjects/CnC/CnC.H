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
    Foam::CnC

Group
    grpFunctionObjects

Description
    This function object executes system calls, entered in the form of a
    string lists.  Calls can be made at the following points in the
    calculation:
    - every time step
    - every output time
    - end of the calculation

    Example of function object specification:
    \verbatim
    systemCall1
    {
        type        systemCall;
        functionObjectLibs ("libsystemCall.so");
        ...
        executeCalls
        (
            "echo execute"
        );
        writeCalls
        (
            "echo \*\*\* writing data \*\*\*"
        );
        endCalls
        (
            "echo \*\*\* writing .bashrc \*\*\*"
            "cat ~/.bashrc"
            "echo \*\*\* done \*\*\*"
        );
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: systemCall   | yes         |
        executeCalls | list of calls on execute | yes        |
        writeCalls   | list of calls on write  | yes         |
        endCalls     | list of calls on end    | yes         |
    \endtable

Note
    Since this function object executes system calls, there is a potential
    security risk.  In order to use the \c systemCall function object, the
    \c allowSystemOperations must be set to '1'; otherwise, system calls will
    not be allowed.

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    systemCall.C
    IOsystemCall.H

\*---------------------------------------------------------------------------*/

#ifndef CnC_H
#define CnC_H

#include "stringList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                         Class CnC Declaration
\*---------------------------------------------------------------------------*/

class CnC
{
protected:

    // Private data

        //- Name of this set of system calls
        word name_;
		
		//- Reference to the database
		const objectRegistry& obr_;

        //- List of calls to execute - every step
        stringList executeCalls_;

        //- List of calls to execute when exiting the time-loop
        stringList endCalls_;

        //- List of calls to execute - write steps
        stringList writeCalls_;
		
		//- List of calls to execute - write steps (only when restart is requested)
		stringList restartCalls_;
		
		bool restartRequested_;
		
		scalar maxClockTime_;
		scalar deltaT_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CnC(const CnC&);

        //- Disallow default bitwise assignment
        void operator=(const CnC&);


public:

    //- Runtime type information
    TypeName("CnC");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        CnC
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary&,
            const bool loadFromFilesUnused = false
        );


    //- Destructor
    virtual ~CnC();


    // Member Functions

        //- Return name of the system call set
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the system calls
        virtual void read(const dictionary&);

        //- Execute the "executeCalls" at each time-step
        virtual void execute();

        //- Execute the "endCalls" at the final time-loop
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write, execute the "writeCalls"
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
