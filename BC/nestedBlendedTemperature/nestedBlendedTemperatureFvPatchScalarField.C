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

#include "nestedBlendedTemperatureFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nestedBlendedTemperatureFvPatchScalarField::
nestedBlendedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{}


Foam::nestedBlendedTemperatureFvPatchScalarField::
nestedBlendedTemperatureFvPatchScalarField
(
    const nestedBlendedTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::nestedBlendedTemperatureFvPatchScalarField::
nestedBlendedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );
}


Foam::nestedBlendedTemperatureFvPatchScalarField::
nestedBlendedTemperatureFvPatchScalarField
(
    const nestedBlendedTemperatureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf)
{}


Foam::nestedBlendedTemperatureFvPatchScalarField::
nestedBlendedTemperatureFvPatchScalarField
(
    const nestedBlendedTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::nestedBlendedTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    const label& patchIndex(this->patch().index());
    const volScalarField& srcField =
            this->db().lookupObject<volScalarField>("TDsBlended");
    const scalarField& srcPatch(srcField.boundaryField()[patchIndex]);
    this->operator==(srcPatch);
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::nestedBlendedTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       nestedBlendedTemperatureFvPatchScalarField
   );
}


// ************************************************************************* //
