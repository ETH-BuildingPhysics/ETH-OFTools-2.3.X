/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
   \\/     M anipulation  |
-------------------------------------------------------------------------------
License
  This file is part of OpenFOAM.
  OpenFOAM is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 3 of the License, or (at your
  option) any later version.
  OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.
  You should have received a copy of the GNU General Public License
  along with OpenFOAM; if not, write to the Free Software Foundation,
  Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
\*---------------------------------------------------------------------------*/

#include "inflowGenerator.H"
#include "addToRunTimeSelectionTable.H"
#include "Vector.H"
#include "SubField.H"
#include "IFstream.H"
#include "OFstream.H"
#include "boundBox.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inflowGenerator::
inflowGenerator
(
  const fvPatch& p,
  const DimensionedField<vector, volMesh>& iF
)
:
  fixedValueFvPatchField<vector>(p, iF),
  LY_(0.0),
  LZ_(0.0),
  dy_(0.0),
  dz_(0.0),
  NY_(0),
  NZ_(0),
  gridFactor_(1.0),
  curTimeIndex_(-1),
  UMean(p.size()),
  ReStress(p.size()),
  uFluctFiltered(p.size()),
  uFluctTemporal_old(p.size()),
  uFluctTemporal(p.size()),
  uFluctFinal(p.size()),
  origin_(vector::zero),
  mapperVP_Ptr_(NULL),
  mapperIV_Ptr_(NULL),
  mapperIP_Ptr_(NULL),
  perturb_(0),
  Lund_(),
  ranGen(label(time(0))),
  cleanRestart_(false),
  isInitialized_(false),
  isRestart_(false),
  correlationShape_("exp")
{}


Foam::inflowGenerator::
inflowGenerator
(
  const inflowGenerator& ptf,
  const fvPatch& p,
  const DimensionedField<vector, volMesh>& iF,
  const fvPatchFieldMapper& mapper
)
:
  fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
  LY_(ptf.LY_),
  LZ_(ptf.LZ_),
  dy_(ptf.dy_),
  dz_(ptf.dz_),
  NY_(ptf.NY_),
  NZ_(ptf.NZ_),
  gridFactor_(ptf.gridFactor_),
  curTimeIndex_(-1),
  UMean(ptf.UMean,mapper),
  ReStress(ptf.ReStress,mapper),
  uFluctFiltered(ptf.uFluctFiltered,mapper),
  uFluctTemporal_old(ptf.uFluctTemporal_old,mapper),
  uFluctTemporal(ptf.uFluctTemporal,mapper),
  uFluctFinal(ptf.uFluctFinal,mapper),
  origin_(ptf.origin_),
  mapperVP_Ptr_(NULL),
  mapperIV_Ptr_(NULL),
  mapperIP_Ptr_(NULL),
  perturb_(ptf.perturb_),
  Lund_(ptf.Lund_,mapper),
  ranGen(label(time(0))),
  cleanRestart_(ptf.cleanRestart_),
  isInitialized_(ptf.isInitialized_),
  isRestart_(ptf.isRestart_),
  correlationShape_(ptf.correlationShape_)
{}


Foam::inflowGenerator::
inflowGenerator
(
  const fvPatch& p,
  const DimensionedField<vector, volMesh>& iF,
  const dictionary& dict
)
:
  fixedValueFvPatchField<vector>(p, iF, dict),
  gridFactor_(dict.lookupOrDefault("gridFactor", 1.0)),
  curTimeIndex_(-1),
  UMean(p.size(),pTraits<vector>::zero),
  ReStress(p.size(),pTraits<symmTensor>::zero),
  uFluctFiltered(p.size(),pTraits<vector>::zero),
  uFluctTemporal_old(p.size(),pTraits<vector>::zero),
  uFluctTemporal(p.size(),pTraits<vector>::zero),
  uFluctFinal(p.size(),pTraits<vector>::zero),
  //origin_(dict.lookup("origin")),
  origin_(vector::zero),
  mapperVP_Ptr_(NULL),
  mapperIV_Ptr_(NULL),
  mapperIP_Ptr_(NULL),
  perturb_(dict.lookupOrDefault("perturb", 1e-6)),
  Lund_(p.size(), pTraits<tensor>::zero),
  ranGen(label(time(0))),
  cleanRestart_(false),
  isInitialized_(false),
  isRestart_(false),
  correlationShape_(dict.lookupOrDefault<word>("correlationShape", "exp"))
{
  //Info << "Foam::inflowGenerator dict constructor" << endl;

  if (dict.found("cleanRestart"))
  {
    cleanRestart_ = true;
  }



  if (dict.found("uFluctTemporal") && !cleanRestart_)
  {
    isRestart_=true;
    uFluctTemporal = vectorField("uFluctTemporal", dict, p.size());
  }
  Info << "Correlation: "<< correlationShape_ << endl;
  //Info << "constructor finished" << endl;
}


Foam::inflowGenerator::
inflowGenerator
(
  const inflowGenerator& ptf
)
:
  fixedValueFvPatchField<vector>(ptf),
  LY_(ptf.LY_),
  LZ_(ptf.LZ_),
  dy_(ptf.dy_),
  dz_(ptf.dz_),
  NY_(ptf.NY_),
  NZ_(ptf.NZ_),
  gridFactor_(ptf.gridFactor_),
  curTimeIndex_(-1),
  UMean(ptf.UMean),
  ReStress(ptf.ReStress),
  uFluctFiltered(ptf.uFluctFiltered),
  uFluctTemporal_old(ptf.uFluctTemporal_old),
  uFluctTemporal(ptf.uFluctTemporal),
  uFluctFinal(ptf.uFluctFinal),
  origin_(ptf.origin_),
  mapperVP_Ptr_(NULL),
  mapperIV_Ptr_(NULL),
  mapperIP_Ptr_(NULL),
  perturb_(ptf.perturb_),
  Lund_(ptf.Lund_),
  ranGen(ptf.ranGen),
  cleanRestart_(ptf.cleanRestart_),
  isInitialized_(ptf.isInitialized_),
  isRestart_(ptf.isRestart_),
  correlationShape_(ptf.correlationShape_)
{}


Foam::inflowGenerator::
inflowGenerator
(
  const inflowGenerator& ptf,
  const DimensionedField<vector, volMesh>& iF
)
:
  fixedValueFvPatchField<vector>(ptf, iF),
  LY_(ptf.LY_),
  LZ_(ptf.LZ_),
  dy_(ptf.dy_),
  dz_(ptf.dz_),
  NY_(ptf.NY_),
  NZ_(ptf.NZ_),
  gridFactor_(ptf.gridFactor_),
  curTimeIndex_(-1),
  UMean(ptf.UMean),
  ReStress(ptf.ReStress),
  uFluctFiltered(ptf.uFluctFiltered),
  uFluctTemporal_old(ptf.uFluctTemporal_old),
  uFluctTemporal(ptf.uFluctTemporal),
  uFluctFinal(ptf.uFluctFinal),
  origin_(ptf.origin_),
  mapperVP_Ptr_(NULL),
  mapperIV_Ptr_(NULL),
  mapperIP_Ptr_(NULL),
  perturb_(ptf.perturb_),
  Lund_(ptf.Lund_),
  ranGen(ptf.ranGen),
  cleanRestart_(ptf.cleanRestart_),
  isInitialized_(ptf.isInitialized_),
  isRestart_(ptf.isRestart_),
  correlationShape_(ptf.correlationShape_)
{}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
/*-------------------------------------------------------------------------------*\
                            FilteredNoiseInflowGenerator
                            using Selective Filtering
\*-------------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * autoSizeGrid  * * * * * * * * * * * * * * //
void Foam::inflowGenerator::autoSizeGrid()
{
  /*
  gets the bounding box of the inlet patch cell centers, and computes the virtual grid
  spacing using the sqrt of the smallest cell size.

  origin_     : min. coordinates of all inlet faces
  dy_ = dz_   : virtual grid spacing
  NY_         : number of grid points in y direction
  NZ_         : number of grid points in z direction
  */

  boundBox localBb(patch().Cf());
  //Info << localBb << endl;
  origin_ = localBb.min();
  LZ_=localBb.span().component(2);
  LY_=localBb.span().component(1);

  scalar smallestFace=0;
  smallestFace=gMin(patch().magSf());
  //Info << "cell dy: "<< Foam::sqrt(smallestFace) << endl;
  scalar dx=Foam::sqrt(smallestFace)*gridFactor_;

  //same grid spacing for y and z direction
  NY_=LY_/dx+1;
  NZ_=LZ_/dx+1;

  dy_=dx;
  dz_=dx;
  Info << "Virtual grid with origin: " << origin_ << " and nr. of points Y:"<< NY_<<" Z: " << NZ_ << endl;
  Info << "Y coordinate:"<< origin_.component(vector::Y)+NY_*dy_<<" Z coordinate: " << origin_.component(vector::Z)+NZ_*dz_ << endl;
}

// * * * * * * * * * * * * * * * initData  * * * * * * * * * * * * * //
void Foam::inflowGenerator::initData()
{
  //set initial seed for rand() for each proc
  srand((Pstream::myProcNo()+1)*time(NULL));

  //automatically size virtual grid
  autoSizeGrid();

  Info << "Virtual grid spacing dy: " << dy_ << " dz: " << dz_ << endl;

  //############################# Create Virtual Grid ##################################
  // Find Number of Faces on inlet patch
  labelList nFacesProc_(Pstream::nProcs(),0);
  nFacesProc_[Pstream::myProcNo()] = this->patch().Cf().size();
  reduce(nFacesProc_,sumOp<labelList>());
  label nFaces_ = sum(nFacesProc_);
  vectorField fCentres_;
  fCentres_.setSize(nFaces_,vector::zero);
  // Get all cellCentres on inlet patch
  List<vectorField> fCentresProc_;
  fCentresProc_.setSize(Pstream::nProcs());
  label sizeF = nFacesProc_[Pstream::myProcNo()];

  label startF = 0;
  for(int i=0;i<Pstream::myProcNo();i++){
    startF = startF + nFacesProc_[i];
  }

  fCentresProc_[Pstream::myProcNo()]=SubField<vector>(fCentres_,sizeF,startF);
  fCentresProc_[Pstream::myProcNo()]=patch().Cf();

  Pstream::gatherList(fCentresProc_);
  Pstream::scatterList(fCentresProc_);

  fCentres_ = ListListOps::combine<vectorField>(fCentresProc_,accessOp<vectorField>());

  if(debug)
  {
    if(Pstream::master())
    {
      fileName rootPath(this->db().time().constant()/"boundaryData"/this->patch().name());
      OFstream(rootPath/"dbg_fCentres_")() << fCentres_;
    }
  }

  //###################################################################################
  // New Sort Method to sort for z and y Coordinates

  // First Find all unique y and z Coordinates
  scalarField zSort, ySort, zSortFinal, ySortFinal;
  zSort.setSize(nFaces_,0);
  ySort.setSize(nFaces_,0);
  scalar zMin, zMax, yMin, yMax;
  zMin = origin_.component(vector::Z)-dz_;
  yMin = origin_.component(vector::Y)-dy_;
  label zCount, yCount;
  zCount=0;
  yCount=0;
  zMax=origin_.component(vector::Z);
  yMax=origin_.component(vector::Y);


  boundBox localBb(patch().Cf());
  vector originOpo_ = localBb.max();

  while(zMax<roundSix(originOpo_.component(vector::Z)))
  {
    zMax = roundSix(originOpo_.component(vector::Z));
    forAll(fCentres_,I)
    {
      if(roundSix(fCentres_[I].component(2))<zMax)
      {
        if(roundSix(fCentres_[I].component(2))>zMin)
        {
          zMax = roundSix(fCentres_[I].component(2));
        }
      }
    }
    zSort[zCount] = zMax;
    zMin = zMax;
    zCount++;
  }
  zSortFinal.setSize(zCount,0);
  forAll(zSortFinal,I)
  {
    zSortFinal[I]=zSort[I];
  }

  while(yMax<roundSix(originOpo_.component(vector::Y)))
  {
    yMax = roundSix(originOpo_.component(vector::Y));
    forAll(fCentres_,I)
    {
      if(roundSix(fCentres_[I].component(1))<yMax)
      {
        if(roundSix(fCentres_[I].component(1))>yMin)
        {
          yMax = roundSix(fCentres_[I].component(1));
        }
      }
    }
    ySort[yCount] = yMax;
    yMin = yMax;
    yCount++;
  }
  ySortFinal.setSize(yCount,0);
  forAll(ySortFinal,I)
  {
    ySortFinal[I]=ySort[I];
  }

  // Now Sort for y and z
  label sortCount=0;
  vectorField fSorted_;
  fSorted_.setSize(nFaces_,vector::zero);
  forAll(zSortFinal,Iz)
  {
    forAll(ySortFinal,Iy)
    {
      forAll(fCentres_,I)
      {
        if(roundSix(fCentres_[I].component(2))==zSortFinal[Iz])
        {
          if(roundSix(fCentres_[I].component(1))==ySortFinal[Iy])
          {
            fSorted_[sortCount].component(0) = roundSix(fCentres_[I].component(0));
            fSorted_[sortCount].component(1) = roundSix(fCentres_[I].component(1));
            fSorted_[sortCount].component(2) = roundSix(fCentres_[I].component(2));
            sortCount++;
          }
        }
      }
    }
  }

  //###################################################################################
  // Delete fCentres at the boundary Layer / bottom of inlet from list of indices
  vectorField fReducedFinal_;
  bool deleteSmallFaces = false;
  if(deleteSmallFaces)
  {
    vectorField fReduced_;
    fReduced_.setSize(nFaces_,vector::zero);
    label rCount=1;
    fReduced_[0]=fSorted_[0];
    label Lcount=0;
    label lineNumber = 2; // After this many horizontal lines of coordinates no further horizontal lines will be skipped...
    for(int i=1;i<nFaces_;i++)
    {
      if(fSorted_[i].component(2)==fReduced_[rCount-1].component(2))
      {
        if(fSorted_[i].component(1)>=(fReduced_[rCount-1].component(1)+dy_*0.95))
        {
          fReduced_[rCount]=fSorted_[i];
          rCount++;
        }
      }
      else if((fSorted_[i].component(2)>=(fReduced_[rCount-1].component(2)+dz_*0.95))||(Lcount>=lineNumber))
      {
        fReduced_[rCount]=fSorted_[i];
        rCount++;
        Lcount++;
      }
    }
    fReducedFinal_.setSize(rCount,vector::zero);
    forAll(fReducedFinal_,I)
    {
      fReducedFinal_[I]=fReduced_[I];
    }
  }
  else
  {
    label rCount = nFaces_;
    fReducedFinal_.setSize(rCount,vector::zero);
    forAll(fReducedFinal_,I)
    {
      fReducedFinal_[I]=fSorted_[I];
    }
  }

  //###################################################################################
  // resort for y and then z for distributed load on Procs:
  // Works if integral length scale are increasing with respect to z
  vectorField fCentresFinal_;
  fCentresFinal_.setSize(fReducedFinal_.size(),vector::zero);
  sortCount = 0;
  label StartSort = 0;
  label SortIndex = 0;
  while(sortCount<fReducedFinal_.size()){
    if(SortIndex<fReducedFinal_.size())
    {
      fCentresFinal_[sortCount] = fReducedFinal_[SortIndex];
      SortIndex += Pstream::nProcs();
      sortCount++;
    }
    else
    {
      StartSort++;
      SortIndex = StartSort;
    }
  }


  //###################################################################################
  // Find nearest virtualGridpoints to cellCentres
  virualGridPoints_.setSize(fCentresFinal_.size(),vector::zero);
  forAll(fCentresFinal_,fI)
  {
    scalar VGPY = 0;
    scalar VGPZ = 0;
    scalar distY = fCentresFinal_[fI].component(1)-origin_.component(vector::Y);
    scalar distZ = fCentresFinal_[fI].component(2)-origin_.component(vector::Z);
    scalar remY = distY-floor(distY/dy_)*dy_;
    scalar remZ = distZ-floor(distZ/dz_)*dz_;
    if (remY>dy_/2) {
      VGPY = origin_.component(vector::Y)+distY-remY+dy_;
    }
    else {
      VGPY = origin_.component(vector::Y)+distY-remY;
    }
    if (remZ>dz_/2) {
      VGPZ = origin_.component(vector::Z)+distZ-remZ+dz_;
    }
    else {
      VGPZ = origin_.component(vector::Z)+distZ-remZ;
    }
    virualGridPoints_[fI] = vector(origin_.component(vector::X),VGPY,VGPZ);
  }

  //###################################################################################


  //########################### Indices for VirtualGrid ################################
  virtualFilteredField_.setSize(fCentresFinal_.size(),vector::zero);
  yindices_.setSize(virtualFilteredField_.size(),0);
  zindices_.setSize(virtualFilteredField_.size(),0);
  forAll(virtualFilteredField_,I)
  {
    yindices_[I] = round((virualGridPoints_[I].component(1)-origin_.component(vector::Y))/dy_);
    zindices_[I] = round((virualGridPoints_[I].component(2)-origin_.component(vector::Z))/dz_);
  }

  //###################################################################################

  if(debug) //debug
  {
    // Find Size of patch Part on each Proc
    forAll(nFacesProc_,I)
    {
      Info << "Proc # " << I << " has PatchSize = " << nFacesProc_[I] << endl;
    }
    if(Pstream::master())
    {
      fileName rootPath(this->db().time().constant()/"boundaryData"/this->patch().name());
      OFstream(rootPath/"dbg_zSortFinal")() << zSortFinal;
      OFstream(rootPath/"dbg_ySortFinal")() << ySortFinal;
      OFstream(rootPath/"dbg_fCentres_")() << fCentres_;
      OFstream(rootPath/"dbg_fSorted_")() << fSorted_;
      OFstream(rootPath/"dbg_fReducedFinal_")() << fReducedFinal_;
      OFstream(rootPath/"dbg_fCentresFinal_")() << fCentresFinal_;
      OFstream(rootPath/"dbg_virualGridPoints_")() << virualGridPoints_;
      OFstream(rootPath/"dbg_yindices_")() << yindices_;
      OFstream(rootPath/"dbg_zindices_")() << zindices_;
    }
  }


  //determine how many virtual grid points are assigned to one processor. If nr. virtual points is not a multiple of nr. processors,
  //the last couple processors may have one less indices than others.
  indicesPerProc_ = 0;
  if(Pstream::master())
  {
    Info << "Generating Inflow for " << virtualFilteredField_.size() << " Faces out of "
    << fCentres_.size() << " ( " << scalar(virtualFilteredField_.size())/fCentres_.size()*100 << "% )" << endl;
    indicesPerProc_ = floor(scalar(virtualFilteredField_.size())/(Pstream::nProcs()));
    rest_=virtualFilteredField_.size()-(indicesPerProc_*(Pstream::nProcs()));
    Info << "Distributing Indices per Proc: " << indicesPerProc_ << endl;
    if(rest_>0)
    {
      Info << "First " << rest_ << " Procs will do +1 Indices" << endl;
    }
  }
  Pstream::scatter(indicesPerProc_);
  Pstream::scatter(rest_);

  //Info << "Foam::inflowGenerator read input data" << endl;

  pointField samplePointsField;
  vectorField refField_;
  Field<symmTensor> RField_;
  vectorField LxFullField_;
  vectorField LyFullField_;
  vectorField LzFullField_;

  if(Pstream::master())
  {
    pointIOField samplePoints
    (
      IOobject
      (
        "points",
        this->db().time().constant(),
        "boundaryData"/this->patch().name(),
        this->db(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        false
      )
    );

    // Read values from File
    vectorIOField ref_
    (
      IOobject
      (
        "ref",
        this->db().time().constant(),
        "boundaryData"
        /this->patch().name(),
        this->db(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        false
      )
    );

    IOField<symmTensor> R_
    (
      IOobject
      (
        "R",
        this->db().time().constant(),
        "boundaryData"
        /this->patch().name(),
        this->db(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        false
      )
    );

    vectorIOField Lx_
    (
      IOobject
      (
        "Lx",
        this->db().time().constant(),
        "boundaryData"
        /this->patch().name(),
        this->db(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        false
      )
    );

    vectorIOField Ly_
    (
      IOobject
      (
        "Ly",
        this->db().time().constant(),
        "boundaryData"
        /this->patch().name(),
        this->db(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        false
      )
    );

    vectorIOField Lz_
    (
      IOobject
      (
        "Lz",
        this->db().time().constant(),
        "boundaryData"
        /this->patch().name(),
        this->db(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        false
      )
    );

    samplePointsField=samplePoints;
    refField_=ref_;
    RField_=R_;
    LxFullField_=Lx_;
    LyFullField_=Ly_;
    LzFullField_=Lz_;
  }

  Pstream::scatter(samplePointsField);
  Pstream::scatter(refField_);
  Pstream::scatter(RField_);
  Pstream::scatter(LxFullField_);
  Pstream::scatter(LyFullField_);
  Pstream::scatter(LzFullField_);

  //Info << "Foam::inflowGenerator map inflow data" << endl;
  // Allocate the interpolator
  mapperIP_Ptr_.reset
  (
    new pointToPointPlanarInterpolation
    (
      samplePointsField,
      this->patch().patch().faceCentres(),
      perturb_
    )
  );
  //Info << "Foam::inflowGenerator mapperIP_Ptr_ created" << endl;

  UMean = mapperIP_Ptr_().interpolate(refField_);
  ReStress = mapperIP_Ptr_().interpolate(RField_);
  //Info << "Foam::inflowGenerator interpolated UMean,ReStress" << endl;

  forAll(ReStress, I)
  {
    if ((ReStress[I].xx() < 0) || (ReStress[I].yy() < 0) || (ReStress[I].zz() < 0))
      FatalErrorIn("inflowGenerator::")<< "Some RMS in input file R are negative"<<abort(FatalError);
  }

  //fill the lund tensor array
  Lund_.replace(tensor::XX, sqrt(ReStress.component(symmTensor::XX)));
  Lund_.replace(tensor::YX, ReStress.component(symmTensor::XY)/Lund_.component(tensor::XX));
  Lund_.replace(tensor::ZX, ReStress.component(symmTensor::XZ)/Lund_.component(tensor::XX));
  Lund_.replace(tensor::YY, sqrt(ReStress.component(symmTensor::YY)-sqr(Lund_.component(tensor::YX))));
  Lund_.replace(tensor::ZY, (ReStress.component(symmTensor::YZ) - Lund_.component(tensor::YX)*Lund_.component(tensor::ZX) )/Lund_.component(tensor::YY));
  Lund_.replace(tensor::ZZ, sqrt(ReStress.component(symmTensor::ZZ) - sqr(Lund_.component(tensor::ZX))-sqr(Lund_.component(tensor::ZY))));
  //Info << "Foam::inflowGenerator Lund_ created" << endl;


  mapperVP_Ptr_.reset
  (
    new pointToPointPlanarInterpolation
    (
      virualGridPoints_,
      this->patch().patch().faceCentres(),
      perturb_,
      true // nearestOnly = true;
    )
  );

  //Info << "Foam::inflowGenerator mapperVP_Ptr_ created" << endl;

  //Interpolate from input data to virtual grid

  mapperIV_Ptr_.reset
  (
    new pointToPointPlanarInterpolation
    (
      samplePointsField,
      virualGridPoints_,
      perturb_
    )
  );

  //Info << "Foam::inflowGenerator mapperIV_Ptr_ created" << endl;

  //map Ly, Lz field to virtual grid
  LyField_ = mapperIV_Ptr_().interpolate(LyFullField_);
  LzField_ = mapperIV_Ptr_().interpolate(LzFullField_);

  //Info << "Foam::inflowGenerator LyField_,LzField_ interpolated" << endl;

  //map Lx field to patch
  LxField_ = mapperIP_Ptr_().interpolate(LxFullField_);

  //Info << "Foam::inflowGenerator LxFullField_ interpolated" << endl;

  //convert length scale to virtual grid points
  NLyField_.setSize(LyField_.size());
  NLzField_.setSize(LzField_.size());



  forAll(LyField_, I)
  {
    NLyField_[I]=vector(round(LyField_[I].component(0)/dy_),round(LyField_[I].component(1)/dy_),round(LyField_[I].component(2)/dy_));
    NLzField_[I]=vector(round(LzField_[I].component(0)/dz_),round(LzField_[I].component(1)/dz_),round(LzField_[I].component(2)/dz_));
  }

  //convert x-length scale to timestep-scale (using frozen turbulence assumption)
  NLxField_.setSize(LxField_.size());
  scalar dt =  db().time().deltaT().value();

  forAll(LxField_, I)
  {
    scalar Ux = UMean[I].component(0);
    label nx = 0;
    label ny = 0;
    label nz = 0;
    if (Ux!=0)
    {
      nx = round((LxField_[I].component(0) / Ux) / dt);
      ny = round((LxField_[I].component(1) / Ux) / dt);
      nz = round((LxField_[I].component(2) / Ux) / dt);
    }
    NLxField_[I]=labelVector(nx,ny,nz);
  }

  //For current processor, get start and end index in virtual list.
  label start, size;
  if (Pstream::myProcNo() < rest_)
  {
    start = Pstream::myProcNo()*indicesPerProc_ + Pstream::myProcNo();
    size = indicesPerProc_+1;
  }
  else
  {
    start = Pstream::myProcNo()*indicesPerProc_ + rest_;
    size = indicesPerProc_;
  }
  labelListList procIdx; //dummy list with correct number of indices for each proc. Can be used to loop over indices per proc.
  procIdx.setSize(Pstream::nProcs());
  procIdx[Pstream::myProcNo()].setSize(size,0);

  //initialize with number of processors
  filterCoeff_yz_u_Proc.setSize(Pstream::nProcs());
  filterCoeff_yz_v_Proc.setSize(Pstream::nProcs());
  filterCoeff_yz_w_Proc.setSize(Pstream::nProcs());

  //only fill element of current proc
  filterCoeff_yz_u_Proc[Pstream::myProcNo()].setSize(size);
  filterCoeff_yz_v_Proc[Pstream::myProcNo()].setSize(size);
  filterCoeff_yz_w_Proc[Pstream::myProcNo()].setSize(size);

  nfK_=4;
  //loop through all indices in one proc, and set filter kernel size and fill with 0.0
  forAll(procIdx[Pstream::myProcNo()],subI)
  {
    int I = subI+start; //index in full array
    filterCoeff_yz_u_Proc[Pstream::myProcNo()][subI].setSize((nfK_*NLzField_[I].component(0)+1)*(nfK_*NLyField_[I].component(0)+1),0.0);
    filterCoeff_yz_v_Proc[Pstream::myProcNo()][subI].setSize((nfK_*NLzField_[I].component(1)+1)*(nfK_*NLyField_[I].component(1)+1),0.0);
    filterCoeff_yz_w_Proc[Pstream::myProcNo()][subI].setSize((nfK_*NLzField_[I].component(2)+1)*(nfK_*NLyField_[I].component(2)+1),0.0);
  }

  //to determine size of virtual grid, use the largest filter kernel size
  NLyMax_u = 0;
  NLyMax_v = 0;
  NLyMax_w = 0;

  NLzMax_u = 0;
  NLzMax_v = 0;
  NLzMax_w = 0;

  //find largest length scale (relates to filter kernel size)
  forAll(NLyField_, I)
  {
    if (NLyField_[I].component(0)>NLyMax_u)
      NLyMax_u=NLyField_[I].component(0);
    if (NLyField_[I].component(1)>NLyMax_v)
      NLyMax_v=NLyField_[I].component(1);
    if (NLyField_[I].component(2)>NLyMax_w)
      NLyMax_w=NLyField_[I].component(2);

    if (NLzField_[I].component(0)>NLzMax_u)
      NLzMax_u=NLzField_[I].component(0);
    if (NLzField_[I].component(1)>NLzMax_v)
      NLzMax_v=NLzField_[I].component(1);
    if (NLzField_[I].component(2)>NLzMax_w)
      NLzMax_w=NLzField_[I].component(2);
  }

  //set size of virtual grid. (remark: not NY_+nfK_*NLyMax_u+1 since one point is shared...)
  virtualRandomField_u_.setSize((NY_+nfK_*NLyMax_u)*(NZ_+nfK_*NLzMax_u));
  virtualRandomField_v_.setSize((NY_+nfK_*NLyMax_v)*(NZ_+nfK_*NLzMax_v));
  virtualRandomField_w_.setSize((NY_+nfK_*NLyMax_w)*(NZ_+nfK_*NLzMax_w));
  zOffsetRnd_u_ = nfK_/2*NLzMax_u;
  yOffsetRnd_u_ = nfK_/2*NLyMax_u;
  zOffsetRnd_v_ = nfK_/2*NLzMax_v;
  yOffsetRnd_v_ = nfK_/2*NLyMax_v;
  zOffsetRnd_w_ = nfK_/2*NLzMax_w;
  yOffsetRnd_w_ = nfK_/2*NLyMax_w;
  Info << "virtulRandomField has size (" <<
  virtualRandomField_u_.size() << ", " <<
  virtualRandomField_v_.size() << ", " <<
  virtualRandomField_w_.size() << ")" << nl << endl;

  Info << "calculating filter coefficients" << endl;

  //each processor calculates its filter kernels
  forAll(procIdx[Pstream::myProcNo()],subI)
  {
    int I = subI+start; //location that belongs to this processor in the full field

    get2DFilterCoeff_New(filterCoeff_yz_u_Proc[Pstream::myProcNo()][subI],NLyField_[I].component(0), NLzField_[I].component(0));
    get2DFilterCoeff_New(filterCoeff_yz_v_Proc[Pstream::myProcNo()][subI],NLyField_[I].component(1), NLzField_[I].component(1));
    get2DFilterCoeff_New(filterCoeff_yz_w_Proc[Pstream::myProcNo()][subI],NLyField_[I].component(2), NLzField_[I].component(2));
  }

  //Info << filterCoeff_yz_u_Proc[Pstream::myProcNo()][10];
  //Info << "done calculating filter coefficients" << endl;

  if(debug)
  {
    fileName rootPath(this->db().time().constant()/"boundaryData"/this->patch().name());
    OFstream(rootPath/"dbg_pointsInput")() << samplePointsField;
    OFstream(rootPath/"dbg_pointsVirtual")() << virualGridPoints_;
    OFstream(rootPath/"dbg_pointsPatch")() << this->patch().Cf();

    OFstream(rootPath/"dbg_LxFullField_")() << LxFullField_;
    OFstream(rootPath/"dbg_LyFullField_")() << LyFullField_;
    OFstream(rootPath/"dbg_LzFullField_")() << LzFullField_;
    OFstream(rootPath/"dbg_LxField_")() << LxField_;
    OFstream(rootPath/"dbg_LyField_")() << LyField_;
    OFstream(rootPath/"dbg_LzField_")() << LzField_;
    OFstream(rootPath/"dbg_NLxField_")() << NLxField_;
    OFstream(rootPath/"dbg_NLyField_")() << NLyField_;
    OFstream(rootPath/"dbg_NLzField_")() << NLzField_;

    OFstream(rootPath/"dbg_ReStress")() << ReStress;
    OFstream(rootPath/"dbg_UMean")() << UMean;

    OFstream(rootPath/"dbg_filterCoeff_yz_u")() << filterCoeff_yz_u_Proc;
    OFstream(rootPath/"dbg_filterCoeff_yz_v")() << filterCoeff_yz_v_Proc;
    OFstream(rootPath/"dbg_filterCoeff_yz_w")() << filterCoeff_yz_w_Proc;
  }
}

// * * * * * * * * * * * * * * * roundSix  * * * * * * * * * * * * * //
Foam::scalar Foam::inflowGenerator::roundSix(scalar numberN)
{
  // Helper Function to round to 6 decimals
  scalar factor = 1.0*pow(10,6);
  scalar numberR = round(numberN*factor)/factor;
  return numberR;
}


// * * * * * * * * * * * * * * * getRandomNumber  * * * * * * * * * * * * * //
Foam::scalar Foam::inflowGenerator::getRandomNumber()
{
  // return distribution(generator); // is slower

  // generates approximate normal distribution, values between -6/+6
  scalar val=0;
  for(int i=1; i<=12; ++i)
  {
    val += (scalar(rand()) / (scalar(RAND_MAX)+1.0));
  }
  return val-6.0;
}

// * * * * * * * * * * * * * * * getRandomField  * * * * * * * * * * * * * //
void Foam::inflowGenerator::getRandomField()
{
  // Generate Random Field in parallel
  List<scalarField> virtualRandomField_u_Proc_;
  List<scalarField> virtualRandomField_v_Proc_;
  List<scalarField> virtualRandomField_w_Proc_;
  virtualRandomField_u_Proc_.setSize(Pstream::nProcs());
  virtualRandomField_v_Proc_.setSize(Pstream::nProcs());
  virtualRandomField_w_Proc_.setSize(Pstream::nProcs());
  label sizeRFu = floor(virtualRandomField_u_.size()/Pstream::nProcs());
  label sizeRFv = floor(virtualRandomField_v_.size()/Pstream::nProcs());
  label sizeRFw = floor(virtualRandomField_w_.size()/Pstream::nProcs());
  label restRFu = virtualRandomField_u_.size()-sizeRFu*Pstream::nProcs();
  label restRFv = virtualRandomField_v_.size()-sizeRFv*Pstream::nProcs();
  label restRFw = virtualRandomField_w_.size()-sizeRFw*Pstream::nProcs();
  // Get Start and End Index for current Proc
  label startRFu, startRFv, startRFw;
  if(Pstream::myProcNo()<restRFu)
  {
    startRFu = Pstream::myProcNo()*sizeRFu + Pstream::myProcNo();
    sizeRFu++;
  }
  else
  {
    startRFu = Pstream::myProcNo()*sizeRFu + restRFu;
  }

  if(Pstream::myProcNo()<restRFv)
  {
    startRFv = Pstream::myProcNo()*sizeRFv + Pstream::myProcNo();
    sizeRFv++;
  }
  else
  {
    startRFv = Pstream::myProcNo()*sizeRFv + restRFv;
  }

  if(Pstream::myProcNo()<restRFw)
  {
    startRFw = Pstream::myProcNo()*sizeRFw + Pstream::myProcNo();
    sizeRFw++;
  }
  else
  {
    startRFw = Pstream::myProcNo()*sizeRFw + restRFw;
  }

  // Create SubFields
  virtualRandomField_u_Proc_[Pstream::myProcNo()]=SubField<scalar>(virtualRandomField_u_,sizeRFu,startRFu);
  virtualRandomField_v_Proc_[Pstream::myProcNo()]=SubField<scalar>(virtualRandomField_v_,sizeRFv,startRFv);
  virtualRandomField_w_Proc_[Pstream::myProcNo()]=SubField<scalar>(virtualRandomField_w_,sizeRFw,startRFw);

  //Fill random fields
  forAll(virtualRandomField_u_Proc_[Pstream::myProcNo()], i)
  {
    virtualRandomField_u_Proc_[Pstream::myProcNo()][i] = scalar(getRandomNumber());
  }
  forAll(virtualRandomField_v_Proc_[Pstream::myProcNo()], i)
  {
    virtualRandomField_v_Proc_[Pstream::myProcNo()][i] = scalar(getRandomNumber());
  }
  forAll(virtualRandomField_w_Proc_[Pstream::myProcNo()], i)
  {
    virtualRandomField_w_Proc_[Pstream::myProcNo()][i] = scalar(getRandomNumber());
  }
  // Combine all parts to complete random field
  Pstream::gatherList(virtualRandomField_u_Proc_);
  Pstream::scatterList(virtualRandomField_u_Proc_);
  Pstream::gatherList(virtualRandomField_v_Proc_);
  Pstream::scatterList(virtualRandomField_v_Proc_);
  Pstream::gatherList(virtualRandomField_w_Proc_);
  Pstream::scatterList(virtualRandomField_w_Proc_);
  virtualRandomField_u_ = ListListOps::combine<scalarField>(virtualRandomField_u_Proc_,accessOp<scalarField>());
  virtualRandomField_v_ = ListListOps::combine<scalarField>(virtualRandomField_v_Proc_,accessOp<scalarField>());
  virtualRandomField_w_ = ListListOps::combine<scalarField>(virtualRandomField_w_Proc_,accessOp<scalarField>());
}

// * * * * * * * * * * * * * * * getFilterCoeff_New  * * * * * * * * * * * * * //
//1d filter coeff, new format
void Foam::inflowGenerator::getFilterCoeff_New(scalarList& b_x, label NLX_x)
{
  const scalar pi = constant::mathematical::pi;
  if (NLX_x==0)
  {
    b_x.setSize(1);
    b_x[0]=0.0;
  }
  else
  {
    double sumx = 0.0;
    label NLX2P1_x=nfK_*NLX_x+1;
    b_x.setSize(NLX2P1_x);

    for (int j=0; j<NLX2P1_x; j++)
    {
      if (correlationShape_=="exp")
        sumx += Foam::exp(-2.0*j/NLX_x);
      else if (correlationShape_=="doubleExp")
        sumx += Foam::exp(-4.0*fabs((j-(nfK_/2.0)*NLX_x)/(NLX_x)));
      else if (correlationShape_=="gaussian")
        sumx += Foam::exp(-2.0*pi*Foam::sqr(scalar(j-(nfK_/2.0)*NLX_x))/(2.0*Foam::sqr(scalar(NLX_x))));
      else
        Info << "correlationShape" << correlationShape_ << "does not exist (ERROR)" << endl;
    }
    sumx = Foam::sqrt(sumx);

    for (int j=0; j<NLX2P1_x; j++)
    {
      if (correlationShape_=="exp")
        b_x[j] = Foam::exp(-1.0*j/NLX_x)/sumx;
      else if (correlationShape_=="doubleExp")
        b_x[j] = Foam::exp(-2.0*fabs((j-(nfK_/2.0)*NLX_x)/(NLX_x)))/sumx;
      else if (correlationShape_=="gaussian")
        b_x[j]= Foam::exp(-pi*Foam::sqr(scalar(j-(nfK_/2.0)*NLX_x))/(2.0*Foam::sqr(scalar(NLX_x))))/sumx;
      else
        Info << "correlationShape" << correlationShape_ << "does not exist (ERROR)" << endl;
    }
  }
}

// * * * * * * * * * * * * * * * get2DFilterCoeff_New  * * * * * * * * * * * * * //
//2d filter coeff, new format
void Foam::inflowGenerator::get2DFilterCoeff_New(scalarList& filter,label NLY_x, label NLZ_x)
{
  scalarList by_x;
  scalarList bz_x;

  //1D filter kernel
  getFilterCoeff_New(by_x, NLY_x);
  getFilterCoeff_New(bz_x, NLZ_x);

  //2D filter kernel from 1D kernels
  for (int i=0; i<(nfK_*NLY_x+1); i++)
  {
    for (int j=0; j<(nfK_*NLZ_x+1); j++)
    {
      filter[get1DIndex(i,j,nfK_*NLZ_x+1)] = by_x[i]*bz_x[j];
    }
  }
}

// * * * * * * * * * * * * * * * spatialCorr  * * * * * * * * * * * * * //
void Foam::inflowGenerator::spatialCorr()
{
  label start, size;
  if (Pstream::myProcNo() < rest_)
  {
    start = Pstream::myProcNo()*indicesPerProc_ + Pstream::myProcNo();
    size = indicesPerProc_+1;
  }
  else
  {
    start = Pstream::myProcNo()*indicesPerProc_ + rest_;
    size = indicesPerProc_;
  }


  List<vectorField> virtualFilteredFieldProc_;
  virtualFilteredFieldProc_.setSize(Pstream::nProcs());

  virtualFilteredFieldProc_[Pstream::myProcNo()]=SubField<vector>(virtualFilteredField_,size,start);

  bool doSpatialCorr=true;
  //apply filter
  forAll(virtualFilteredFieldProc_[Pstream::myProcNo()],subI)
  {
    if(doSpatialCorr)
    {
      int I = subI+start;
      int i = yindices_[I]; // i = yindices on virtual Grid
      int j = zindices_[I]; // j = zindices on virtual Grid

      double sumTmp_u=0;
      double sumTmp_v=0;
      double sumTmp_w=0;

      for (int ii=0;ii<(nfK_*NLyField_[I].component(0)+1);ii++)
      {
        int start_rnd=get1DIndex(i+yOffsetRnd_u_-nfK_/2*NLyField_[I].component(0)+ii, j+zOffsetRnd_u_-nfK_/2*NLzField_[I].component(0), (NZ_+nfK_*NLzMax_u));
        int size_rnd=nfK_*NLzField_[I].component(0)+1;
        int size_filt=nfK_*NLzField_[I].component(0)+1;
        int start_filt=get1DIndex(ii, 0, (nfK_*NLzField_[I].component(0)+1));

        SubField<scalar> rnd = SubField<scalar>(virtualRandomField_u_,size_rnd,start_rnd);
        SubField<scalar> filt = SubField<scalar>(filterCoeff_yz_u_Proc[Pstream::myProcNo()][subI],size_filt,start_filt);

        sumTmp_u+=sumProd(rnd,filt);
      }

      for (int ii=0;ii<(nfK_*NLyField_[I].component(1)+1);ii++)
      {
        int start_rnd=get1DIndex(i+yOffsetRnd_v_-nfK_/2*NLyField_[I].component(1)+ii, j+zOffsetRnd_v_-nfK_/2*NLzField_[I].component(1), (NZ_+nfK_*NLzMax_v));
        int size_rnd=nfK_*NLzField_[I].component(1)+1;
        int size_filt=nfK_*NLzField_[I].component(1)+1;
        int start_filt=get1DIndex(ii, 0, (nfK_*NLzField_[I].component(1)+1));

        SubField<scalar> rnd = SubField<scalar>(virtualRandomField_v_,size_rnd,start_rnd);
        SubField<scalar> filt = SubField<scalar>(filterCoeff_yz_v_Proc[Pstream::myProcNo()][subI],size_filt,start_filt);

        sumTmp_v+=sumProd(rnd,filt);
      }

      for (int ii=0;ii<(nfK_*NLyField_[I].component(2)+1);ii++)
      {
        int start_rnd=get1DIndex(i+yOffsetRnd_w_-nfK_/2*NLyField_[I].component(2)+ii, j+zOffsetRnd_w_-nfK_/2*NLzField_[I].component(2), (NZ_+nfK_*NLzMax_w));
        int size_rnd=nfK_*NLzField_[I].component(2)+1;
        int size_filt=nfK_*NLzField_[I].component(2)+1;
        int start_filt=get1DIndex(ii, 0, (nfK_*NLzField_[I].component(2)+1));

        SubField<scalar> rnd = SubField<scalar>(virtualRandomField_w_,size_rnd,start_rnd);
        SubField<scalar> filt = SubField<scalar>(filterCoeff_yz_w_Proc[Pstream::myProcNo()][subI],size_filt,start_filt);

        sumTmp_w+=sumProd(rnd,filt);
      }

      virtualFilteredFieldProc_[Pstream::myProcNo()][subI].component(0)=sumTmp_u;
      virtualFilteredFieldProc_[Pstream::myProcNo()][subI].component(1)=sumTmp_v;
      virtualFilteredFieldProc_[Pstream::myProcNo()][subI].component(2)=sumTmp_w;
    }
    else
    {
      virtualFilteredFieldProc_[Pstream::myProcNo()][subI].component(0)=scalar(getRandomNumber());
      virtualFilteredFieldProc_[Pstream::myProcNo()][subI].component(1)=scalar(getRandomNumber());
      virtualFilteredFieldProc_[Pstream::myProcNo()][subI].component(2)=scalar(getRandomNumber());
    }
  }

  Pstream::gatherList(virtualFilteredFieldProc_);
  Pstream::scatterList(virtualFilteredFieldProc_);

  virtualFilteredField_ = ListListOps::combine<vectorField>(virtualFilteredFieldProc_,accessOp<vectorField>());
}


// * * * * * * * * * * * * * * * temporalCorr  * * * * * * * * * * * * * //
void Foam::inflowGenerator::temporalCorr()
{
  forAll(uFluctTemporal, celli)
  {
    label& NLX_u=NLxField_[celli].component(0);
    label& NLX_v=NLxField_[celli].component(1);
    label& NLX_w=NLxField_[celli].component(2);

    if (NLX_u>0)
      uFluctTemporal[celli].component(vector::X) = uFluctTemporal_old[celli].component(vector::X) * Foam::exp(-1.0/NLX_u) +
      uFluctFiltered[celli].component(vector::X) * Foam::sqrt( 1.0 - Foam::exp(-2.0/NLX_u) );
    else
      uFluctTemporal[celli].component(vector::X) = 0.0;

    if (NLX_v>0)
      uFluctTemporal[celli].component(vector::Y) = uFluctTemporal_old[celli].component(vector::Y) * Foam::exp(-1.0/NLX_v) +
      uFluctFiltered[celli].component(vector::Y) * Foam::sqrt( 1.0 - Foam::exp(-2.0/NLX_v) );
    else
      uFluctTemporal[celli].component(vector::Y)=0.0;

    if (NLX_w>0)
      uFluctTemporal[celli].component(vector::Z) = uFluctTemporal_old[celli].component(vector::Z) * Foam::exp(-1.0/NLX_w) +
      uFluctFiltered[celli].component(vector::Z) * Foam::sqrt( 1.0 - Foam::exp(-2.0/NLX_w) );
    else
      uFluctTemporal[celli].component(vector::Z) = 0.0;
  }
  uFluctTemporal_old = uFluctTemporal;
}


// * * * * * * * * * * * * * * * interpolFluct  * * * * * * * * * * * * * //
void Foam::inflowGenerator::interpolFluct()
{
  uFluctFiltered = mapperVP_Ptr_().interpolate(virtualFilteredField_);
}

// * * * * * * * * * * * * * * * scaleFluct  * * * * * * * * * * * * * //
void Foam::inflowGenerator::scaleFluct()
{
  bool doLundScaling=true;
  // Scaling of fluctuations and storing in field uFluctFinal
  // Lund (1998)
  if (doLundScaling)
  {
    uFluctFinal=Lund_&uFluctTemporal;
  }
  else
  {
    uFluctFinal=uFluctTemporal;
  }
}


// * * * * * * * * * * * * * * * massFlowCorr  * * * * * * * * * * * * * //
void Foam::inflowGenerator::massFlowCorr()
{
  scalarField Sf(this->patch().magSf());
  scalarField phi(uFluctFinal.component(0)*this->patch().magSf());
  scalar dPhi = gSum(phi); // global sum (gSum) goes over all procs
  scalar sumSf = gSum(this->patch().magSf());

  Info << "Delta Vdot: " << dPhi << endl;

  forAll(uFluctFinal, celli)
  {
    uFluctFinal[celli].component(0)=uFluctFinal[celli].component(0)-dPhi/sumSf;
  }
  phi=uFluctFinal.component(0)*this->patch().magSf();
  Info << "Corrected Delta Vdot: " << gSum(phi) << endl;
}

// * * * * * * * * * * * * * * * initialize  * * * * * * * * * * * * * //
void Foam::inflowGenerator::initialize()
{
  Info << "FilteredNoiseInflowGenerator using Selective Filtering" << endl;
  timeBu_ = this->db().time().elapsedCpuTime();
  Info << "Initializing Inflow Generation" << endl;
  if (isRestart_ && !cleanRestart_)
  {
    Info << "inflowGenerator: continuing a simulation" << endl;
    uFluctTemporal_old = uFluctTemporal;
    initData();
  }
  else
  {
    // for the very first time step only
    Info << "inflowGenerator: preparing first time step " << endl;
    initData();
    getRandomField();
    spatialCorr();
    interpolFluct();
    uFluctTemporal_old = uFluctFiltered;
  }
}


// * * * * * * * * * * * * * * * get1DIndex  * * * * * * * * * * * * * //
inline int Foam::inflowGenerator::get1DIndex(int x, int y, int yMax)
{
  //helper function to convert 2d array index into 1d list index. Array is converted in y direction into list.
  return x*yMax + y;
}

// * * * * * * * * * * * * * * * get2DIndex  * * * * * * * * * * * * * //
inline void Foam::inflowGenerator::get2DIndex(int I, int& x, int& y, int yMax)
{
  //helper function to get 2d array index (x,y) from 1d list index.
  x=floor(I/yMax);
  y=(I-x*yMax);
}


// * * * * * * * * * * * * * * * updateCoeffs  * * * * * * * * * * * * * //
void Foam::inflowGenerator::updateCoeffs()
{
  if (updated())
  {
    return;
  }

  //initialize data: read input data files, allocate memory, calculate filter coefficients
  if(!isInitialized_)
  {
    initialize();
    isInitialized_=true;
  }

  if (curTimeIndex_ != this->db().time().timeIndex())
  {
    Info << "Starting Inflow Generation" << endl;
    // Info << "Starting Inflow Generation, time = "<<this->db().time().elapsedClockTime()<<" s"<<endl;
    eTime_ = this->db().time().elapsedCpuTime();

    eTimeSub_ = this->db().time().elapsedCpuTime();
    getRandomField();   //update random field
    eTimeSub_ = this->db().time().elapsedCpuTime()-eTimeSub_;
    Info << "getRandomField Time required = " << eTimeSub_ << " s" << endl;
    eTimeSub_ = this->db().time().elapsedCpuTime();
    spatialCorr();      //filter random field
    eTimeSub_ = this->db().time().elapsedCpuTime()-eTimeSub_;
    Info << "spatial Filtering Time required = " << eTimeSub_ << " s" << endl;
    interpolFluct();    //interpolate from virtual grid to patch
    temporalCorr();     //create new temporally correlated slice
    scaleFluct();       //Apply Lund's transformation to new slice
    massFlowCorr();		//Mass Flow Correction

    eTime_ = this->db().time().elapsedCpuTime()-eTime_;
    eTimetot_ = this->db().time().elapsedCpuTime()-timeBu_;
    timeBu_ = this->db().time().elapsedCpuTime();
    // Info << "Finishing Inflow Generation, time = "<<this->db().time().elapsedClockTime()<<" s"<<endl;
    Info << "Finished Inflow Generation" << endl;
    Info << "Time required for Inlet Generation: " << eTime_ << " s" << endl;
    Info << "timeRatio = " << eTime_/eTimetot_ << nl << endl;

    Field<vector>& patchField = *this;
    forAll(patchField, celli)
    {
      patchField[celli] = UMean[celli] + uFluctFinal[celli]; //create final velocity field
    }

    isRestart_=true; //set flag to enable writing of data, allows to restart a simulation from an old timestep
    curTimeIndex_ = this->db().time().timeIndex();
  }

  fixedValueFvPatchField<vector>::updateCoeffs();
}

// * * * * * * * * * * * * * * * write  * * * * * * * * * * * * * //
void Foam::inflowGenerator::write(Ostream& os) const
{
  fvPatchField<vector>::write(os);
  os.writeKeyword("perturb") << perturb_ << token::END_STATEMENT << nl;
  os.writeKeyword("correlationShape") << correlationShape_ << token::END_STATEMENT << nl;
  os.writeKeyword("gridFactor") << gridFactor_ << token::END_STATEMENT << nl;
  //os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
  //os.writeKeyword("NY") << NY_ << token::END_STATEMENT << nl;
  //os.writeKeyword("NZ") << NZ_ << token::END_STATEMENT << nl;
  //os.writeKeyword("SizeY") << LY_ << token::END_STATEMENT << nl;
  //os.writeKeyword("SizeZ") << LZ_ << token::END_STATEMENT << nl;
  writeEntry("value", os);
  if (isRestart_)
    uFluctTemporal.writeEntry("uFluctTemporal", os);

  // debug: write the interpolated source field if required.
  //if (dict.lookupOrDefault("writeSourceFields",false)==true)

  if(debug)
  {
    fileName rootPath(this->db().time().timePath());

    OFstream(rootPath/"dbg_virtualFilteredField")() << virtualFilteredField_;
    OFstream(rootPath/"dbg_uFluctFiltered")() << uFluctFiltered;
    OFstream(rootPath/"dbg_uFluctFinal")() << uFluctFinal;
  }
}

// * * * * * * * * * * * * * * * autoMap  * * * * * * * * * * * * * //
void Foam::inflowGenerator::autoMap(const fvPatchFieldMapper& m)
{
  Field<vector>::autoMap(m);

  uFluctTemporal.autoMap(m);

  // Clear interpolator
  mapperVP_Ptr_.clear();
  mapperIV_Ptr_.clear();
  mapperIP_Ptr_.clear();
}


// * * * * * * * * * * * * * * * rMap  * * * * * * * * * * * * * //
void Foam::inflowGenerator::rmap(const fvPatchField<vector>& ptf, const labelList& addr)
{
  fixedValueFvPatchField<vector>::rmap(ptf, addr);

  const inflowGenerator& tiptf = refCast<const inflowGenerator>(ptf);

  uFluctTemporal.rmap(tiptf.uFluctTemporal, addr);

  // Clear interpolator
  mapperVP_Ptr_.clear();
  mapperIV_Ptr_.clear();
  mapperIP_Ptr_.clear();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  makePatchTypeField
  (
    fvPatchVectorField,
    inflowGenerator
  );
}


// ************************************************************************* //
