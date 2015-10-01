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

Application
    postAvg

Description
    Creates postAvg

Source files:
    createFields.H

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "fvIOoptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("fieldName");

    #include "setRootCase.H"
    #include "createTime.H"

    word fieldName(args.additionalArgs()[0]);

    instantList timeDirs = timeSelector::select0(runTime, args);
	
	runTime.setTime(timeDirs[0], 0);
	
    #include "createMesh.H"
    
	Info<< "    Create mean field" << nl<<endl;
    volScalarField dummy
    (
     IOobject
     (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
      ),
     mesh
    );
	
	// create the new file name with the _mean extension
    const word EXT_MEAN = "Mean";
    const word meanFieldName=fieldName+EXT_MEAN;
    
    volScalarField mean
    (
     IOobject
     (
        meanFieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
      ),
     dummy
    );
	
	// and now initialise at zero. Much better would be to do this right away. 
    // How can that be done?
    mean*=0;

    // start the counter out 0;
    int nfield=0;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        // give some information
        Info<< "timeI: "<< timeI<< " Time = " << runTime.timeName() << endl;


	// read the header field 'fieldName'
        IOobject fieldHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check field exists
        if (fieldHeader.headerOk())
        {
            // the mesh exists, read the data
            mesh.readUpdate();

            if (fieldHeader.headerClassName() == "volScalarField")
            {
               // the data is volScalar data -> add the field to the mean
               Info<< "    Reading volScalarField " << nfield << " " << fieldName << endl;
               volScalarField field(fieldHeader, mesh);
               mean+=field;
               nfield++;
            }
            else
            {
                FatalError
                    << "Only possible to average volScalarFields "
                    << nl << exit(FatalError);
            }
        }
        else
        {
            Info<< "    No field " << fieldName << endl;
        }

        Info<< endl;
    }

	// devide by the number of added fields
    if(nfield>0){
      Info<< "number of fields added: "<< nfield << endl;
      mean/=nfield;
    }

    Info<< "writing to file "  << endl;
    mean.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

