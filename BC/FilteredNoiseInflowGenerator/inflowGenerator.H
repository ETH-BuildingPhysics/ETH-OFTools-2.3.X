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
  Foam::inflowGenerator
Group
  grpInletBoundaryConditions
Description
Note
  correlationShape (string) can be exp or doubleExp or gaussian
SeeAlso
  Foam::fixedValueFvPatchField
SourceFiles
  inflowGenerator.C
\*---------------------------------------------------------------------------*/

#ifndef inflowGenerator_H
#define inflowGenerator_H

#include "fixedValueFvPatchFields.H"
#include "DataEntry.H"
#include "List.H"
#include "labelVector.H"
#include "pointToPointPlanarInterpolation.H"
#include "mathematicalConstants.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
         Class inflowGenerator Declaration
\*---------------------------------------------------------------------------*/

class inflowGenerator
:
  public fixedValueFvPatchVectorField
{
  // Private data

  scalar LY_;
  scalar LZ_;
  scalar dy_;
  scalar dz_;
  label NY_;
  label NZ_;
  scalar gridFactor_;

  //- Current time index (used for updating)
  label curTimeIndex_;

  // Field declaration and initialization
  vectorField UMean;

  Field<symmTensor> ReStress;
  vectorField uFluctFiltered;     // intermediate spatially correlated (filtered) fluctuation on patch mesh
  vectorField uFluctTemporal_old; // old temporal correlated fluctuation
  vectorField uFluctTemporal;     // new temporal correlated fluctuation
  vectorField uFluctFinal;        // final (scaled) fluctuation

  vectorField cells;
  scalarField ycentre;
  scalarField zcentre;

  pointField virualGridPoints_;

  scalarField virtualRandomField_u_;
  scalarField virtualRandomField_v_;
  scalarField virtualRandomField_w_;

  vectorField virtualFilteredField_;

  List<scalarListList> filterCoeff_yz_u_Proc;
  List<scalarListList> filterCoeff_yz_v_Proc;
  List<scalarListList> filterCoeff_yz_w_Proc;

  //Geometry of virtual Grid
  vector origin_;
  scalar hight_;
  scalar width_;

  //Measure Elapsed Time
  scalar eTime_;
  scalar eTimeSub_;
  scalar timeBu_;
  scalar eTimetot_;

  //mapper from virtual grid to patch
  autoPtr<pointToPointPlanarInterpolation> mapperVP_Ptr_;
  //mapper from inlet data to virtual grid
  autoPtr<pointToPointPlanarInterpolation> mapperIV_Ptr_;
  //mapper from inlet data to patch
  autoPtr<pointToPointPlanarInterpolation> mapperIP_Ptr_;
  scalar perturb_;

  Field<tensor> Lund_;
  Foam::Random ranGen;

  vectorField LxField_;
  Field<labelVector> NLxField_;

  vectorField LyField_;
  vectorField LzField_;
  Field<labelVector> NLyField_;
  Field<labelVector> NLzField_;


  label NLyMax_u;
  label NLyMax_v;
  label NLyMax_w;

  label NLzMax_u;
  label NLzMax_v;
  label NLzMax_w;

  //for parallel processing
  label indicesPerProc_;
  label rest_;
  labelList yindices_;
  labelList zindices_;

  //filter kernel factor (multiple of length scale)
  label nfK_;
  label zOffsetRnd_u_;
  label yOffsetRnd_u_;
  label zOffsetRnd_v_;
  label yOffsetRnd_v_;
  label zOffsetRnd_w_;
  label yOffsetRnd_w_;


  bool cleanRestart_;
  bool isInitialized_;
  bool isRestart_;


  word correlationShape_;

  public:

  //- Runtime type information
  TypeName("filteredNoiseInflowGenerator");


  // Constructors

  //- Construct from patch and internal field
  inflowGenerator
  (
    const fvPatch&,
    const DimensionedField<vector, volMesh>&
  );

  //- Construct from patch, internal field and dictionary
  inflowGenerator
  (
    const fvPatch&,
    const DimensionedField<vector, volMesh>&,
    const dictionary&
  );

  //- Construct by mapping given
  //  inflowGenerator
  //  onto a new patch
  inflowGenerator
  (
    const inflowGenerator&,
    const fvPatch&,
    const DimensionedField<vector, volMesh>&,
    const fvPatchFieldMapper&
  );

  //- Construct as copy
  inflowGenerator
  (
    const inflowGenerator&
  );

  //- Construct and return a clone
  virtual tmp<fvPatchVectorField> clone() const
  {
    return tmp<fvPatchVectorField>
    (
      new inflowGenerator(*this)
    );
  }

  //- Construct as copy setting internal field reference
  inflowGenerator
  (
    const inflowGenerator&,
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
      new inflowGenerator(*this, iF)
    );
  }



  // Member functions

  //read experimental data
  void initData();

  //prepare for time loop
  void initialize();

  void autoSizeGrid();
  //generate a normal distributed random number
  virtual scalar getRandomNumber();

  //returns index of a 2D array in list format
  int get1DIndex(int x, int y, int yMax);
  void get2DIndex(int I, int& x, int& y, int yMax);

  // round to 6 decimals
  virtual scalar roundSix(scalar numberN);

  // new get filter coeff
  void getFilterCoeff_New(scalarList& b_x, label NLX_x);
  void get2DFilterCoeff_New(scalarList& filter,label NLY_x, label NLZ_x);

  //generate normal distributed random field
  void getRandomField();

  void spatialCorr();

  //get temporal correlation
  void temporalCorr();

  //interpolating velocity fluctuation from virtual grid on OpenFOAM grid
  //void interpolFluct_bilinear();
  void interpolFluct();

  //getting scaled velocity fluctuation with method of LUND (1998)
  void scaleFluct();

  //Conduct mass Flow correction
  void massFlowCorr();

  virtual void autoMap(const fvPatchFieldMapper&);

  virtual void rmap(const fvPatchField<vector>&, const labelList&);

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
