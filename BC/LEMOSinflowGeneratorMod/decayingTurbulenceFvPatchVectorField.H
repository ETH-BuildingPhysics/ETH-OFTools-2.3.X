/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) held by original author
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
    Foam::decayingTurbulenceFvPatchVectorField

Description
    direction: It is +1 or −1 depending on wheather the flow goes in positive or negative x-direction,respectively.
    LField: The integral length field.
    refField: The external velocity field.
    RField: The Reynolds stress tensor field.
    value: The value of the field; at the initial time it should be equal to the refField field.

    inlet
    {
	type decayingTurbulenceInflowGenerator;
 	direction 1;
	LField nonuniform List<scalar> ...;
	refField nonuniform List<vector> ...;
	RField nonuniform List<symmTensor> ...;
	value nonuniform List<vector> ...;
    }

SourceFiles
    decayingTurbulenceFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef decayingTurbulenceFvPatchVectorField_H
#define decayingTurbulenceFvPatchVectorField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"
#include "SLList.H"
#include "pointToPointPlanarInterpolation.H"

#include "decayingVorton.H"

namespace Foam
{

  
class decayingTurbulenceFvPatchVectorField : public fixedValueFvPatchField<vector>
{
private:

    Field<scalar>          LField_;
    Field<vector>          refField_;
    Field<symmTensor>      RField_;
    label                  curTimeIndex_;
    SLList<decayingVorton> vortons_;
    Field<tensor>          Lund_;
    Field<symmTensor>      R_;
    label                  ind_;
    label                  direction_;

    //- 2D interpolation
    autoPtr<pointToPointPlanarInterpolation> mapperPtr_;
    //- Fraction of perturbation (fraction of bounding box) to add
    scalar perturb_;

public:

    TypeName("decayingTurbulenceInflowGeneratorMod");

    decayingTurbulenceFvPatchVectorField
    (
        const fvPatch&, 
        const DimensionedField<vector, volMesh>&
    );

    decayingTurbulenceFvPatchVectorField
    (
        const fvPatch&, 
        const DimensionedField<vector, volMesh>&, 
        const dictionary&
    );

    decayingTurbulenceFvPatchVectorField
    (
        const decayingTurbulenceFvPatchVectorField&, 
        const fvPatch&, 
        const DimensionedField<vector, volMesh>&, 
        const fvPatchFieldMapper&
    );

    decayingTurbulenceFvPatchVectorField
    (
        const decayingTurbulenceFvPatchVectorField&
    );

    decayingTurbulenceFvPatchVectorField
    (
        const decayingTurbulenceFvPatchVectorField&, 
        const DimensionedField<vector, volMesh>&
    );
    
    virtual tmp<fvPatchField<vector> > clone() const
    {
        return tmp<fvPatchField<vector> >(new decayingTurbulenceFvPatchVectorField(*this));
    }
    
    virtual tmp<fvPatchField<vector> > clone(const DimensionedField<vector, volMesh>& iF) const
    {
        return tmp<fvPatchField<vector> >(new decayingTurbulenceFvPatchVectorField(*this, iF));
    }

    virtual void autoMap(const fvPatchFieldMapper&);

    virtual void rmap(const fvPatchField<vector>&, const labelList&);

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
    
private:
    
    void doUpdate();
    void readInletData();
    void writeSourceFields();

};
} 

#endif
