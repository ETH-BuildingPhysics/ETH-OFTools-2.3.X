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
#include "FieldFunctions.H"
#include "boundBox.H"

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
    isRestart_(false)
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
    isRestart_(ptf.isRestart_)
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
    perturb_(dict.lookupOrDefault("perturb", 1e-5)),
    Lund_(p.size(), pTraits<tensor>::zero),
    ranGen(label(time(0))),
    cleanRestart_(false),
    isInitialized_(false),
    isRestart_(false)
{
    //Info << "Foam::inflowGenerator dict constructor" << endl;

    if (dict.found("cleanRestart"))
    {
        cleanRestart_ = true;
    }

    //number of points on virtual grid

/*
    NZ_ = readLabel(dict.lookup("NZ"));
    NY_ = readLabel(dict.lookup("NY"));

    LZ_ = readScalar(dict.lookup("SizeZ"));
    LY_ = readScalar(dict.lookup("SizeY"));
*/


    if (dict.found("uFluctTemporal") && !cleanRestart_)
    {
        isRestart_=true;
        uFluctTemporal = vectorField("uFluctTemporal", dict, p.size());

    }

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
    isRestart_(ptf.isRestart_)
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
    isRestart_(ptf.isRestart_)
{}

void Foam::inflowGenerator::autoSizeGrid()
{

    boundBox localBb(patch().Cf());
    //Info << localBb << endl;
    origin_ = localBb.min();
    LZ_=localBb.span().component(2);
    LY_=localBb.span().component(1);

    scalar smallestFace=0;
    smallestFace=gMin(patch().magSf());
    //Info << "cell dy: "<< Foam::sqrt(smallestFace) << endl;
    scalar dx=Foam::sqrt(smallestFace);
    origin_.component(1)=origin_.component(1);
    origin_.component(2)=origin_.component(2);
    //LY_=LY_;
    //LZ_=LZ_;

    NY_=LY_/dx+1;
    NZ_=LZ_/dx+1;

    dy_=dx;
    dz_=dx;
    Info << "Virtual grid with origin: " << origin_ << " and nr. of points Y:"<< NY_<<" Z: " << NZ_ << endl;
    Info << "Y coordinate:"<< origin_.component(vector::Y)+NY_*dy_<<" Z coordinate: " << origin_.component(vector::Z)+NZ_*dz_ << endl;
    /*
    int numsf = 0;

    forAll(sf, I)
    {
        if (sf[I] > LField_[I]*LField_[I])
        numsf++;
    }

    if (numsf > 0)
        Pout<<"Warning: there are "<<numsf<<" inlet cells (out of "<<patch().size()<<") where the integral length is smaller than the cell size."<<endl;*/
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::inflowGenerator::initData()
{
    //set initial seed for rand() for each proc
    srand((Pstream::myProcNo()+1)*time(NULL));

    //automatically size virtual grid
    autoSizeGrid();

    //NY * NZ
    //dy_ = LY_/NY_; // grid size in the y-dir of the virtual uniform mesh
    //dz_ = LZ_/NZ_; // grid size in the z-dir of the virtual uniform mesh

    Info << "Virtual grid spacing dy: " << dy_ << " dz: " << dz_ << endl;

    //perturb_=1e-5;//dy_/2;

    //create coordinates of virtual grid
    virualGridPoints_.setSize(NY_*NZ_);
    for(int z=0;z<NZ_;z++)
    {
        for(int y=0;y<NY_;y++)
        {
            virualGridPoints_[get1DIndex(y, z, NZ_)] = vector(origin_.component(vector::X),origin_.component(vector::Y)+y*dy_,origin_.component(vector::Z)+z*dz_);
        }
    }
    //Info<< virualGridPoints_ << endl;

    virtualFilteredField_.setSize(NY_*NZ_,vector::zero);

    //mapping of indices for parallel processing

    yindices_.setSize(virtualFilteredField_.size(),0);
    zindices_.setSize(virtualFilteredField_.size(),0);
    forAll(virtualFilteredField_,I)
    {
        get2DIndex(I, yindices_[I], zindices_[I], NZ_);
    }

    indicesPerProc_ = 0;
    if(Pstream::master())
    {
        Info << "Total Size: " << virtualFilteredField_.size() << endl;
        indicesPerProc_ = floor(scalar(virtualFilteredField_.size())/(Pstream::nProcs()));
        label rest =0;
        rest=virtualFilteredField_.size()-(indicesPerProc_*(Pstream::nProcs()));
        //Info << indicesPerProc << " rest: " << rest << endl;
        if(rest>0)
        {
            indicesPerProc_ = floor(scalar(virtualFilteredField_.size())/(Pstream::nProcs()-1));
            rest=virtualFilteredField_.size()-(indicesPerProc_*(Pstream::nProcs()-1));
        }

        Info << indicesPerProc_ << " rest: " << rest << endl;

    }
    Pstream::scatter(indicesPerProc_);

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
            perturb_
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

    //convert length scale to grid points
    NLyField_.setSize(LyField_.size());
    NLzField_.setSize(LzField_.size());

    NLxField_.setSize(LxField_.size());

    forAll(LyField_, I)
    {
        NLyField_[I]=vector(round(LyField_[I].component(0)/dy_),round(LyField_[I].component(1)/dy_),round(LyField_[I].component(2)/dy_));
        NLzField_[I]=vector(round(LzField_[I].component(0)/dz_),round(LzField_[I].component(1)/dz_),round(LzField_[I].component(2)/dz_));
    }

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


    label start=Pstream::myProcNo()*indicesPerProc_;
    label size=indicesPerProc_;
    if ((start+size)>virtualFilteredField_.size())
    {
        size=virtualFilteredField_.size()-start;
    }

    //Pout << "start: " << start << " size: " << size << endl;

    filterCoeff_yz_u_Proc.setSize(Pstream::nProcs());
    filterCoeff_yz_v_Proc.setSize(Pstream::nProcs());
    filterCoeff_yz_w_Proc.setSize(Pstream::nProcs());

    labelListList procIdx;
    procIdx.setSize(Pstream::nProcs());
    procIdx[Pstream::myProcNo()].setSize(size,0);

    filterCoeff_yz_u_Proc[Pstream::myProcNo()].setSize(size);
    filterCoeff_yz_v_Proc[Pstream::myProcNo()].setSize(size);
    filterCoeff_yz_w_Proc[Pstream::myProcNo()].setSize(size);

    forAll(procIdx[Pstream::myProcNo()],subI)
    {
        int I = subI+start;
        filterCoeff_yz_u_Proc[Pstream::myProcNo()][subI].setSize((2*NLzField_[I].component(0)+1)*(2*NLyField_[I].component(0)+1),0.0);
        filterCoeff_yz_v_Proc[Pstream::myProcNo()][subI].setSize((2*NLzField_[I].component(1)+1)*(2*NLyField_[I].component(1)+1),0.0);
        filterCoeff_yz_w_Proc[Pstream::myProcNo()][subI].setSize((2*NLzField_[I].component(2)+1)*(2*NLyField_[I].component(2)+1),0.0);
    }

    //set size of filters, get max size
    //max sizes

    NLyMax_u = 0;
    NLyMax_v = 0;
    NLyMax_w = 0;

    NLzMax_u = 0;
    NLzMax_v = 0;
    NLzMax_w = 0;

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

    //Info << NLyMax_u << maximum << endl;

    //set size of virtual grid
    virtualRandomField_u_.setSize((NY_+2*NLyMax_u)*(NZ_+2*NLzMax_u));
    virtualRandomField_v_.setSize((NY_+2*NLyMax_v)*(NZ_+2*NLzMax_v));
    virtualRandomField_w_.setSize((NY_+2*NLyMax_w)*(NZ_+2*NLzMax_w));

    Info << "calculating filter coefficients" << endl;

    forAll(procIdx[Pstream::myProcNo()],subI)
    {
        int I = subI+start;

        get2DFilterCoeff_New(filterCoeff_yz_u_Proc[Pstream::myProcNo()][subI],NLyField_[I].component(0), NLzField_[I].component(0));
        get2DFilterCoeff_New(filterCoeff_yz_v_Proc[Pstream::myProcNo()][subI],NLyField_[I].component(1), NLzField_[I].component(1));
        get2DFilterCoeff_New(filterCoeff_yz_w_Proc[Pstream::myProcNo()][subI],NLyField_[I].component(2), NLzField_[I].component(2));
    }
	
    Info << filterCoeff_yz_u_Proc[Pstream::myProcNo()][10];
    //Info << "done calculating filter coefficients" << endl;
}


Foam::scalar Foam::inflowGenerator::getRandomNumber()
{
    // generates approximate normal distribution, values between -6/+6
    scalar val=0;
    for(int i=1; i<=12; ++i)
    {
        val += (scalar(rand()) / (scalar(RAND_MAX)+1.0));
    }
    return val-6.0;
}


void Foam::inflowGenerator::getRandomField()
{
    if(Pstream::master())
    {
        //Fill random fields
        forAll(virtualRandomField_u_, i)
        {
                virtualRandomField_u_[i] = scalar(getRandomNumber());
        }
        forAll(virtualRandomField_v_, i)
        {
                virtualRandomField_v_[i] = scalar(getRandomNumber());
        }
        forAll(virtualRandomField_w_, i)
        {
                virtualRandomField_w_[i] = scalar(getRandomNumber());
        }
    }
    Pstream::scatter(virtualRandomField_u_);
    Pstream::scatter(virtualRandomField_v_);
    Pstream::scatter(virtualRandomField_w_);
}

//1d filter coeff, new format
void Foam::inflowGenerator::getFilterCoeff_New(scalarList& b_x, label NLX_x)
{
    double sumx = 0.0;
    label NLX2P1_x=2*NLX_x+1;

    b_x.setSize(NLX2P1_x);

    for (int j=0; j<NLX2P1_x; j++)
    {
        //sumx += Foam::sqr(Foam::exp(-2.0*fabs((j-NLX_x)/(NLX_x))));
        sumx += Foam::exp(-4.0*j/NLX_x);
    }
    sumx = Foam::sqrt(sumx);

    for (int j=0; j<NLX2P1_x; j++)
    {
        //b_x[j] = Foam::exp(-2.0*fabs((j-NLX_x)/(NLX_x)))/sumx;
        b_x[j] = Foam::exp(-2.0*j/NLX_x)/sumx;
    }
}

//1d filter coeff, new format
void Foam::inflowGenerator::get2DFilterCoeff_New(scalarList& filter,label NLY_x, label NLZ_x)
{
    scalarList by_x;
    scalarList bz_x;

    getFilterCoeff_New(by_x, NLY_x);
    getFilterCoeff_New(bz_x, NLZ_x);

    for (int i=0; i<(2*NLY_x+1); i++)
    {
        for (int j=0; j<(2*NLZ_x+1); j++)
        {
            filter[get1DIndex(i,j,2*NLZ_x+1)] = by_x[i]*bz_x[j];
        }
    }
}


void Foam::inflowGenerator::spatialCorr()
{

    label start=Pstream::myProcNo()*indicesPerProc_;
    label size=indicesPerProc_;
    if ((start+size)>virtualFilteredField_.size())
    {
        size=virtualFilteredField_.size()-start;
    }

    List<vectorField> virtualFilteredFieldProc_;
    virtualFilteredFieldProc_.setSize(Pstream::nProcs());

    virtualFilteredFieldProc_[Pstream::myProcNo()]=SubField<vector>(virtualFilteredField_,size,start);

    bool doSpatialCorr=true;
    //appy filter
    forAll(virtualFilteredFieldProc_[Pstream::myProcNo()],subI)
    {
        if(doSpatialCorr)
        {
            int I = subI+start;
            int i = yindices_[I];
            int j = zindices_[I];

            double sumTmp_u=0;
            double sumTmp_v=0;
            double sumTmp_w=0;

            for (int ii=0;ii<(2*NLyField_[I].component(0)+1);ii++)
            {
                int start_rnd=get1DIndex(i+NLyMax_u-NLyField_[I].component(0)+ii, j+NLzMax_u-NLzField_[I].component(0), (NZ_+2*NLzMax_u));
                int size_rnd=2*NLzField_[I].component(0)+1;
                int size_filt=2*NLzField_[I].component(0)+1;
                int start_filt=get1DIndex(ii, 0, (2*NLzField_[I].component(0)+1));

                SubField<scalar> rnd = SubField<scalar>(virtualRandomField_u_,size_rnd,start_rnd);
                SubField<scalar> filt = SubField<scalar>(filterCoeff_yz_u_Proc[Pstream::myProcNo()][subI],size_filt,start_filt);

                sumTmp_u+=sumProd(rnd,filt);
                //sumTmp_u+=rnd[0];
            }

            for (int ii=0;ii<(2*NLyField_[I].component(1)+1);ii++)
            {
                int start_rnd=get1DIndex(i+NLyMax_v-NLyField_[I].component(1)+ii, j+NLzMax_v-NLzField_[I].component(1), (NZ_+2*NLzMax_v));
                int size_rnd=2*NLzField_[I].component(1)+1;
                int size_filt=2*NLzField_[I].component(1)+1;
                int start_filt=get1DIndex(ii, 0, (2*NLzField_[I].component(1)+1));

                SubField<scalar> rnd = SubField<scalar>(virtualRandomField_v_,size_rnd,start_rnd);
                SubField<scalar> filt = SubField<scalar>(filterCoeff_yz_v_Proc[Pstream::myProcNo()][subI],size_filt,start_filt);

                sumTmp_v+=sumProd(rnd,filt);
                //sumTmp_v=virtualRandomField_v_[start_rnd];
            }

            for (int ii=0;ii<(2*NLyField_[I].component(2)+1);ii++)
            {
                int start_rnd=get1DIndex(i+NLyMax_w-NLyField_[I].component(2)+ii, j+NLzMax_w-NLzField_[I].component(2), (NZ_+2*NLzMax_v));
                int size_rnd=2*NLzField_[I].component(2)+1;
                int size_filt=2*NLzField_[I].component(2)+1;
                int start_filt=get1DIndex(ii, 0, (2*NLzField_[I].component(2)+1));

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

/*
void Foam::inflowGenerator::interpolFluct_bilinear()
{
    forAll(uFluctFiltered, celli)
    {
        double y_coord = 0.0;
        double z_coord = 0.0;
        int jj = 0;
        int kk = 0;

        for(int j=0; j<NY; j++)
        {
            if (dy + origin_.y() + j*dy > ycentre[celli])
            {
                y_coord = dy + origin_.y() + (j-1)*dy;
                jj=j-1;
                break;
            }
        }


        for(int k=0; k<NZ; k++)
        {
            if (dz + origin_.z() + k*dz > zcentre[celli])
            {
                z_coord = dz + origin_.z() + (k-1)*dz;
                kk=k-1;
                break;
            }
        }

        double y1 = y_coord;
        double y2 = y_coord + dy;
        double z1 = z_coord;
        double z2 = z_coord + dz;

        if (y1==ycentre[celli] && z1==zcentre[celli])
        {
            uFluctFiltered[celli].component(vector::X) = uFluctSpatial[jj][kk];
            uFluctFiltered[celli].component(vector::Y) = vFluctSpatial[jj][kk];
            uFluctFiltered[celli].component(vector::Z) = wFluctSpatial[jj][kk];
        }
        else if (y2==ycentre[celli] && z2==zcentre[celli])
        {
            uFluctFiltered[celli].component(vector::X) = uFluctSpatial[jj+1][kk+1];
            uFluctFiltered[celli].component(vector::Y) = vFluctSpatial[jj+1][kk+1];
            uFluctFiltered[celli].component(vector::Z) = wFluctSpatial[jj+1][kk+1];
        }
        else if (y1==ycentre[celli] && z2==zcentre[celli])
        {
            uFluctFiltered[celli].component(vector::X) = uFluctSpatial[jj][kk+1];
            uFluctFiltered[celli].component(vector::Y) = vFluctSpatial[jj][kk+1];
            uFluctFiltered[celli].component(vector::Z) = wFluctSpatial[jj][kk+1];
        }
        else if (y2==ycentre[celli] && z1==zcentre[celli])
        {
            uFluctFiltered[celli].component(vector::X) = uFluctSpatial[jj+1][kk];
            uFluctFiltered[celli].component(vector::Y) = vFluctSpatial[jj+1][kk];
            uFluctFiltered[celli].component(vector::Z) = wFluctSpatial[jj+1][kk];
        }
        else
        {
            // Bi-linear interpolation
            // for every OpenFoam face patch:
                  // get center of face patch
                  // find 4 neighbouring cell points on virtual grid (y_coord & z_coord)
                  // get fluctuation values of 4 neighbouring cells
                  // first linear interpolation to find fluctuatin of R1 and R2
                  // second linear interpolation to find final fluctuation
            double R1[3] = {0.0, 0.0, 0.0};
            double R2[3] = {0.0, 0.0, 0.0};

            R1[0] = (y2-ycentre[celli])/(y2-y1)*uFluctSpatial[jj][kk] + (ycentre[celli]-y1)/(y2-y1)*uFluctSpatial[jj+1][kk];
            R2[0] = (y2-ycentre[celli])/(y2-y1)*uFluctSpatial[jj][kk+1] + (ycentre[celli]-y1)/(y2-y1)*uFluctSpatial[jj+1][kk+1];
            R1[1] = (y2-ycentre[celli])/(y2-y1)*vFluctSpatial[jj][kk] + (ycentre[celli]-y1)/(y2-y1)*vFluctSpatial[jj+1][kk];
            R2[1] = (y2-ycentre[celli])/(y2-y1)*vFluctSpatial[jj][kk+1] + (ycentre[celli]-y1)/(y2-y1)*vFluctSpatial[jj+1][kk+1];
            R1[2] = (y2-ycentre[celli])/(y2-y1)*wFluctSpatial[jj][kk] + (ycentre[celli]-y1)/(y2-y1)*wFluctSpatial[jj+1][kk];
            R2[2] = (y2-ycentre[celli])/(y2-y1)*wFluctSpatial[jj][kk+1] + (ycentre[celli]-y1)/(y2-y1)*wFluctSpatial[jj+1][kk+1];

            uFluctFiltered[celli].component(vector::X) = (z2-zcentre[celli])/(z2-z1)*R1[0]
                    + (zcentre[celli]-z1)/(z2-z1)*R2[0];
            uFluctFiltered[celli].component(vector::Y) = (z2-zcentre[celli])/(z2-z1)*R1[1]
                    + (zcentre[celli]-z1)/(z2-z1)*R2[1];
            uFluctFiltered[celli].component(vector::Z) = (z2-zcentre[celli])/(z2-z1)*R1[2]
                    + (zcentre[celli]-z1)/(z2-z1)*R2[2];
        }
    }
}
*/

void Foam::inflowGenerator::interpolFluct()
{
    uFluctFiltered = mapperVP_Ptr_().interpolate(virtualFilteredField_);
}


void Foam::inflowGenerator::scaleFluct()
{
    bool doLundScaling=true;
    // Scaling of fluctuations and storing in field uFluctFinal
    // Lund (1998)
    if (doLundScaling)
    {
        uFluctFinal=uFluctTemporal&Lund_;
    }else
    {
        uFluctFinal=uFluctTemporal;
    }
}


void Foam::inflowGenerator::initialize()
{	
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

inline int Foam::inflowGenerator::get1DIndex(int x, int y, int yMax)
{
    //helper function to convert 2d array index into 1d list index
    return x*yMax + y;
}

inline void Foam::inflowGenerator::get2DIndex(int I, int& x, int& y, int yMax)
{
    //helper function to get 2d index (x,y) from 1d index
    x=floor(I/yMax);
    y=(I-x*yMax);
}



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
        Info << "Starting Inflow Generation, time = "<<this->db().time().elapsedClockTime()<<" s"<<endl;


        getRandomField();   //update random field
        spatialCorr();      //filter random field
        interpolFluct();    //interpolate from virtual grid to patch
        temporalCorr();     //create new temporally correlated slice
        scaleFluct();       //Apply Lund's transformation to new slice

        Info << "Finishing Inflow Generation, time = "<<this->db().time().elapsedClockTime()<<" s"<<endl;

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


void Foam::inflowGenerator::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("NY") << NY_ << token::END_STATEMENT << nl;
    os.writeKeyword("NZ") << NZ_ << token::END_STATEMENT << nl;
    os.writeKeyword("SizeY") << LY_ << token::END_STATEMENT << nl;
    os.writeKeyword("SizeZ") << LZ_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
    if (isRestart_)
            uFluctTemporal.writeEntry("uFluctTemporal", os);
}

void Foam::inflowGenerator::autoMap(const fvPatchFieldMapper& m)
{
    Field<vector>::autoMap(m);

    uFluctTemporal.autoMap(m);

    // Clear interpolator
    mapperVP_Ptr_.clear();
    mapperIV_Ptr_.clear();
    mapperIP_Ptr_.clear();

}

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
