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

#include "nestedBlendedVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nestedBlendedVelocityFvPatchVectorField::
nestedBlendedVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{}


Foam::nestedBlendedVelocityFvPatchVectorField::
nestedBlendedVelocityFvPatchVectorField
(
    const nestedBlendedVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{}


Foam::nestedBlendedVelocityFvPatchVectorField::
nestedBlendedVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF)
{
    fvPatchField<vector>::operator=
    (
        vectorField("value", dict, p.size())
    );
}


Foam::nestedBlendedVelocityFvPatchVectorField::
nestedBlendedVelocityFvPatchVectorField
(
    const nestedBlendedVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf)
{}


Foam::nestedBlendedVelocityFvPatchVectorField::
nestedBlendedVelocityFvPatchVectorField
(
    const nestedBlendedVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::nestedBlendedVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

//    const word& patchName(this->patch().name());
    const label& patchIndex(this->patch().index());

    const volVectorField& srcField =
            this->db().lookupObject<volVectorField>("UDsBlended");

    const vectorField& srcPatch(srcField.boundaryField()[patchIndex]);

    this->operator==(srcPatch);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::nestedBlendedVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       nestedBlendedVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
