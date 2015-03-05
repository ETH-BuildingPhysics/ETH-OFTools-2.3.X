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
    UName_("U")
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

    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField* QPtr
        (
            new volScalarField
            (
                IOobject
                (
                    type(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("0", dimless/sqr(dimTime), 0.0)
            )
        );

        mesh.objectRegistry::store(QPtr);
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

        const volVectorField& U =
            mesh.lookupObject<volVectorField>(UName_);

        const volTensorField gradU(fvc::grad(U));

        volScalarField& Q =
            const_cast<volScalarField&>
            (
                mesh.lookupObject<volScalarField>(type())
            );

        Q = 0.5*(sqr(tr(gradU)) - tr(((gradU) & (gradU))));
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
        const volScalarField& Q =
            obr_.lookupObject<volScalarField>(type());

        Info<< type() << " " << name_ << " output:" << nl
            << "    writing field " << Q.name() << nl
            << endl;

        Q.write();
    }
}


// ************************************************************************* //
