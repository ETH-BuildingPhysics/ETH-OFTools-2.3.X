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

Class
    Foam::nestedBlendedVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Velocity boundary condition which should be used with the solver
    nestedDsPimpleFoam. This BC should be used with every boundary listed under
    "blendedPatches" in "netingDict" of nestedDsPimpleFoam

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    nestedBlendedVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef nestedBlendedVelocityFvPatchVectorField_H
#define nestedBlendedVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
           Class nestedBlendedVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class nestedBlendedVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data


public:

   //- Runtime type information
   TypeName("nestedBlendedVelocity");


   // Constructors

        //- Construct from patch and internal field
        nestedBlendedVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nestedBlendedVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nestedBlendedVelocityFvPatchVectorField
        //  onto a new patch
        nestedBlendedVelocityFvPatchVectorField
        (
            const nestedBlendedVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nestedBlendedVelocityFvPatchVectorField
        (
            const nestedBlendedVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new nestedBlendedVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nestedBlendedVelocityFvPatchVectorField
        (
            const nestedBlendedVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new nestedBlendedVelocityFvPatchVectorField(*this, iF)
            );
        }

    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
