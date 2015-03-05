/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "addScalarTimeU.H"
#include "volFields.H"
#include "dictionary.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(addScalarTimeU, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::addScalarTimeU::addScalarTimeU
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    UName_("U"),
    scalarFields_(dict.lookup("scalarFields"))
//    scalarFields_()
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "addScalarTimeU::addScalarTimeU"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

    read(dict);
    Info << "*** scalarFields_ = " << scalarFields_ << endl;

    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        const volScalarField& Sa = mesh.lookupObject<volScalarField>("Sa");
        const dimensionSet& SaDim = Sa.dimensions();
        Info << "*** SaDim = " << SaDim << endl;

        const volVectorField& U = mesh.lookupObject<volVectorField>(UName_);
        const dimensionSet& UDim = U.dimensions();
        Info << "*** UDim = " << UDim << endl;

        volVectorField* STimeUPtr
        (
            new volVectorField
            (
                IOobject
                (
                    "SaU",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("SaU", SaDim*UDim, vector::zero)
            )
        );

        mesh.objectRegistry::store(STimeUPtr);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::addScalarTimeU::~addScalarTimeU()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::addScalarTimeU::read(const dictionary& dict)
{
    if (active_)
    {
        UName_ = dict.lookupOrDefault<word>("UName", "U");
    }
}


void Foam::addScalarTimeU::execute()
{
    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        const volVectorField& U = mesh.lookupObject<volVectorField>(UName_);

        const volScalarField& S = mesh.lookupObject<volScalarField>("Sa");

//        const volTensorField gradU(fvc::grad(U));

//        volScalarField& Q =
//            const_cast<volScalarField&>
//            (
//                mesh.lookupObject<volScalarField>(type())
//            );

//        Q = 0.5*(sqr(tr(gradU)) - tr(((gradU) & (gradU))));

        volVectorField& STimeU =
            const_cast<volVectorField&>
            (
                mesh.lookupObject<volVectorField>("SaU")
            );

        STimeU = S*U;
    }
}


void Foam::addScalarTimeU::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::addScalarTimeU::timeSet()
{
    // Do nothing
}


void Foam::addScalarTimeU::write()
{
    if (active_)
    {
        const volVectorField& STimeU =
            obr_.lookupObject<volVectorField>("SaU");

        Info<< type() << " " << name_ << " output:" << nl
            << "    writing field " << STimeU.name() << nl
            << endl;

        STimeU.write();
    }
}


// ************************************************************************* //
