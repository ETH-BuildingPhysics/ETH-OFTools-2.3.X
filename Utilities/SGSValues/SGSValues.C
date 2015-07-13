/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    cellVolumes

Description
    Calculates and writes the cell volumes  mesh.V() (scalar)

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    IOobject nuSgsHeader
    (
        "nuSgs",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    
    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (Uheader.headerOk())
    {
        Info<< "    Reading U" << endl;
        volVectorField U(Uheader, mesh);
        volTensorField gradU(fvc::grad(U));
        volScalarField nuSgs(nuSgsHeader,mesh);
    

        //volTensorField gradU("gradU",fvc::grad(U)());
    
        volScalarField G("kSgsProd", 2.0*nuSgs*magSqr(symm(gradU)));

        Info<< "Writing G to " << G.name() << " in "
         << runTime.timeName() << endl;

        G.write();
        
}
    Info<< "End" << endl;
}


// ************************************************************************* //
