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

#include "scalarVelocityProduct.H"
#include "volFields.H"
#include "dictionary.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(scalarVelocityProduct, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarVelocityProduct::scalarVelocityProduct
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
    scalarFields_(dict.lookup("scalarFields")),
    USNames_(scalarFields_.size()),
    writeFields_(false)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "scalarVelocityProduct::scalarVelocityProduct"
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

        const volVectorField& U = mesh.lookupObject<volVectorField>(UName_);
        const dimensionSet& UDim = U.dimensions();

        forAll (scalarFields_, i)
        {
            const volScalarField& S =
                mesh.lookupObject<volScalarField>(scalarFields_[i]);
            const dimensionSet& SDim = S.dimensions();

            volVectorField* USproductPtr
            (
                new volVectorField
                (
                    IOobject
                    (
                        USNames_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedVector
                    (
                        USNames_[i]+"Dim", UDim*SDim, vector::zero
                    )
                )
            );

            mesh.objectRegistry::store(USproductPtr);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::scalarVelocityProduct::~scalarVelocityProduct()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::scalarVelocityProduct::read(const dictionary& dict)
{
    if (active_)
    {
        UName_ = dict.lookupOrDefault<word>("UName", "U");
        writeFields_ = dict.lookupOrDefault<bool>("writeFields", false);
        forAll (scalarFields_, i)
        {
            USNames_[i] = "U"+scalarFields_[i];
        }
    }
}


void Foam::scalarVelocityProduct::execute()
{
    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        const volVectorField& U = mesh.lookupObject<volVectorField>(UName_);

        forAll (scalarFields_, i)
        {
            const volScalarField& S =
                mesh.lookupObject<volScalarField>(scalarFields_[i]);

            volVectorField& USproduct =
                const_cast<volVectorField&>
                (
                    mesh.lookupObject<volVectorField>(USNames_[i])
                );

            USproduct = U*S;
        }
    }
}


void Foam::scalarVelocityProduct::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::scalarVelocityProduct::timeSet()
{
    // Do nothing
}


void Foam::scalarVelocityProduct::write()
{
    if (active_ && writeFields_)
    {
        Info<< type() << " " << name_ << " output:" << endl;
        forAll (scalarFields_, i)
        {
            const volVectorField& USproduct =
                obr_.lookupObject<volVectorField>(USNames_[i]);
            Info<< "    writing field " << USproduct.name() << endl;

            USproduct.write();
        }
        Info << nl << endl;
    }
}


// ************************************************************************* //
