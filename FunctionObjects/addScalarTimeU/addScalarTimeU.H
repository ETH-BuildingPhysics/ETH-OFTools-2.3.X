/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::addScalarTimeU

Group
    grpUtilitiesFunctionObjects

Description
    For each given scalar field Si, this function object stores the Si*U fields
    on the mesh database for further manipulation. After averaging such fields,
    The fields SiPrimeUPrimeMean can be calculated with the ad hoc
    post-processing utility.

    Work Flow example:
        During the run:
            1. Calculate the scalar transport of Sa and Sb with the function
               object scalarTransport
            2. Use this function object to add SaU and SbU to the mesh database
               and write the if necessary.
            3. Average SaU and SbU with the function object fieldAverage.
        After the run:
            1. Use the post-processing tool "toolName" to get SaPrimeUPrimeMean
               and SbPrimeUPrimeMean

    addScalarTimeU
    {
        type                 addScalarTimeU;
        functionObjectLibs   ("libaddScalarTimeU.so");
        enabled              true;
        writeSUfields        false;
        outputControl        outputTime;
        region               region0;
        scalarFields
        (
            Sa
            T
            k
        );
    }

SourceFiles
    addScalarTimeU.C

\*---------------------------------------------------------------------------*/

#ifndef addScalarTimeU_H
#define addScalarTimeU_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "OFstream.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                          Class addScalarTimeU Declaration
\*---------------------------------------------------------------------------*/

class addScalarTimeU
{
    // Private data

        //- Name of this set of addScalarTimeU objects
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Name of velocity field, default is "U"
        word UName_;

        //- list of scalar field
        wordList scalarFields_;

        //- list of names of the scalarTimeU fields
        wordList SUNames_;

        //- bool to control writing. default=false
        bool writeSUfields_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        addScalarTimeU(const addScalarTimeU&);

        //- Disallow default bitwise assignment
        void operator=(const addScalarTimeU&);


public:

    //- Runtime type information
    TypeName("addScalarTimeU");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        addScalarTimeU
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~addScalarTimeU();


    // Member Functions

        //- Return name of the set of addScalarTimeU
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the addScalarTimeU data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the Q and write
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
