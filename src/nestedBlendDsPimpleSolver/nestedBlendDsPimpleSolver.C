/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "nestedBlendDsPimpleSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
//        defineTypeNameAndDebug(nestedBlendDsPimpleSolver, 0);
//        defineRunTimeSelectionTable(nestedBlendDsPimpleSolver, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nestedBlendDsPimpleSolver::nestedBlendDsPimpleSolver
(
    word name,
    Time& runTime,
    bool runNested
)
:
    name_(name),
    runTime_(runTime),
    mesh_
    (
        Foam::IOobject
        (
            name_,
            runTime_.timeName(),
            runTime_,
            Foam::IOobject::MUST_READ
        )
    ),
    nestingControlDict_
    (
        IOobject
        (
            "nestingControlDict",
            runTime_.system(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    p_
    (
      IOobject
      (
          "p",
          runTime_.timeName(),
          mesh_,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ),
      mesh_
    ),
    U_
    (
      IOobject
      (
          "U",
          runTime_.timeName(),
          mesh_,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ),
      mesh_
    ),
    phi_
    (
        IOobject
        (
            "phi",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U_) & mesh_.Sf()
    ),
    functionObjects_(runTime_,nestingControlDict_,true),
    laminarTransport_(U_, phi_),
    turbulence_(incompressible::turbulenceModel::New(U_, phi_, laminarTransport_)),
    fvOptions_(mesh_),
    pimple_(mesh_),
    pRefCell_(0),
    pRefValue_(0.0),
    cumulativeContErr_(0),
    nestedCases_(nestingControlDict_.lookup("nestedCases")),
    nestedSolvers_(),
    runNested_(runNested),
    nameMst_(nestingControlDict_.lookup("masterName")),
    meshMst_(runTime_.lookupObject<fvMesh>(nameMst_)),
    interpolator_(meshMst_,mesh_,nestingControlDict_),
    UMstOnNtd_
    (
        IOobject
        (
            "UMstOnNtd",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("velocity",dimLength/dimTime,vector::zero)
    ),
    U0MstOnNtd_
    (
        IOobject
        (
            "U0MstOnNtd",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("velocity",dimLength/dimTime,vector::zero)
    ),
    UDsBlended_
    (
        IOobject
        (
            "UDsBlended",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("velocity",dimLength/dimTime,vector::zero)
    ),
    pMstOnNtd_
    (
        IOobject
        (
            "pMstOnNtd",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("pressure",p_.dimensions(),0)
    ),
    p0MstOnNtd_
    (
        IOobject
        (
            "p0MstOnNtd",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("pressure",p_.dimensions(),0)
    ),
    pDsBlended_
    (
        IOobject
        (
            "pDsBlended",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("pressure",p_.dimensions(),0)
    ),
    wDs_
    (
        IOobject
        (
            "wDs",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("scalarZero",dimless,0),
        "zeroGradient"                         // zeroGradient at the boundary
    )

{
    Info << "Create mesh for time = " << runTime_.timeName() << nl << endl;
    Info << "Reading field p\n" << endl;
    Info << "Reading field U\n" << endl;
    Info << "Reading/calculating face flux field phi\n" << endl;

    setRefCell(p_, mesh_.solutionDict().subDict("PIMPLE"), pRefCell_, pRefValue_);
    interpolator_.creatInterpolator();
    setwDs();
    functionObjects_.start();

    // write wDs in the inital timestep folder.
    wDs_.write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

nestedBlendDsPimpleSolver::~nestedBlendDsPimpleSolver()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void nestedBlendDsPimpleSolver::setwDs()
{
    const wordReList blendedPatches(nestingControlDict_.lookup("blendedPatches"));
    const dimensionedScalar blendingDist
    (
        "lenght",
        dimLength,
        nestingControlDict_.lookup("blendingDist")
    );

    if (blendedPatches.size()>0)
    {
        labelHashSet blendedPatchSet
        (
            mesh_.boundaryMesh().patchSet(blendedPatches)
        );

        patchWave wave
        (
            mesh_,
            blendedPatchSet,
            true
        );

        scalarField& wDsIn = wDs_.internalField();
        wDsIn = wave.distance();

        // linear blending
        forAll (wDsIn, ccI)
        {
            wDsIn[ccI] =
                       - (1.0/blendingDist.value())
                       * wDsIn[ccI] + 1.0;
            if (wDsIn[ccI] < 0.0)
            {
                wDsIn[ccI] = 0.0;
            }
        }

        // set exactly blendingStrength add the faceCell
        forAll (blendedPatchSet, patchID)
        {
            forAll (mesh_.boundary()[patchID], facei)
            {
                wDs_[mesh_.boundary()[patchID].faceCells()[facei]] = 1.0;
            }
        }
    }

    // scale wDs and wUs
    scalar DsBlendingStrength
    (
        nestingControlDict_.lookupOrDefault<scalar>("DsBlendingStrength",1.0)
    );
    wDs_ = DsBlendingStrength*wDs_;

    // apply zeroGradient after internalField modifications
    wDs_.correctBoundaryConditions();
}


void nestedBlendDsPimpleSolver::interpolateMstFields()
{
    // interpolate master U on UMstOnNtd and master U_0 on U0MstOnNtd.
    // As the interpolation procedure involved important interprocessor
    // communication, U0MstOnNtd is simply copied from UMstOnNtd before
    // Updating this last one. On startup, U0MstOnNtd must be interpolatated
    // as UMstOnNtd doesn't still exist (if loop).
    if (U_.timeIndex()-runTime_.time().startTimeIndex() == 1)
    {
        interpolator_.DSinterpolatorObj().mapSrcToTgt
        (
            meshMst_.lookupObject<volVectorField>("U").oldTime(),
            eqOp<vector>(),
            U0MstOnNtd_
        );
        interpolator_.DSinterpolatorObj().mapSrcToTgt
        (
            meshMst_.lookupObject<volScalarField>("p").oldTime(),
            eqOp<scalar>(),
            p0MstOnNtd_
        );
    }
    else
    {
        U0MstOnNtd_ = UMstOnNtd_;
        p0MstOnNtd_ = pMstOnNtd_;
    }

    interpolator_.DSinterpolatorObj().mapSrcToTgt
    (
        meshMst_.lookupObject<volVectorField>("U"),
        eqOp<vector>(),
        UMstOnNtd_
    );

    interpolator_.DSinterpolatorObj().mapSrcToTgt
    (
        meshMst_.lookupObject<volScalarField>("p"),
        eqOp<scalar>(),
        pMstOnNtd_
    );
}


void nestedBlendDsPimpleSolver::setNestedSolvers
(
    List< nestedBlendDsPimpleSolver* > nestedSolvers
)
{
    nestedSolvers_.clear();
    forAll (nestedSolvers,i)
    {
        forAll(nestedCases_,wi)
        {
            if (nestedSolvers[i]->name()==nestedCases_[wi])
            {
                nestedSolvers_.append(nestedSolvers[i]);
            }
        };
    };
}


void nestedBlendDsPimpleSolver::runFunctionObjects()
{
    dimensionedScalar timeValue("timeValue", dimTime, runTime_.value());
    bool running = timeValue < (runTime_.endTime() - 0.1*runTime_.deltaT());

    if (!running)
    {
        // Note, end() also calls an indirect start() as required
        functionObjects_.end();
    }
    else
    {
        functionObjects_.execute();
    }

    running = timeValue < (runTime_.endTime() - 0.1*runTime_.deltaT());
}


void nestedBlendDsPimpleSolver::solve()
{
    dimensionedScalar DsRelaxTime
    (
        "time",
        dimTime,
        nestingControlDict_.lookup("DsRelaxationTime")
    );
    const label& subDeltaT = readLabel(nestingControlDict_.lookup("subDeltaT"));

    // Enter in subCycle time loop for the nested simulation
    scalar subCycleIndex(0.0);
    interpolateMstFields();

    for
    (
        subCycleTime ntdSubCycle
        (
            const_cast<Time&>(runTime_),
            subDeltaT
        );
        !(++ntdSubCycle).end();
    )
    {
        scalar subTimeFrac = (subCycleIndex+1)/subDeltaT;
        UDsBlended_ = (1-subTimeFrac)*U0MstOnNtd_ + subTimeFrac*UMstOnNtd_;
        pDsBlended_ = (1-subTimeFrac)*p0MstOnNtd_ + subTimeFrac*pMstOnNtd_;

        Info << name_ << ": Time = " << runTime_.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple_.loop())
        {
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple_.correct())
            {
                #include "pEqn.H"
            }

            if (pimple_.turbCorr())
            {
                turbulence_->correct();
            }
        }

        #include "CourantNo.H"

        Info << "ExecutionTime = " << runTime_.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime_.elapsedClockTime() << " s"
             << nl << endl;

        runFunctionObjects();

        if (nestedCases_.size()!=0 && runNested_==true)
        {
            forAll (nestedSolvers_, i)
            {
                nestedSolvers_[i]->solve();
            }
        }

        subCycleIndex = subCycleIndex+1.0;
    }
}

// * * * * * * * * * * * * * * * Access Functions  * * * * * * * * * * * * * //
word& nestedBlendDsPimpleSolver::name()
{
    return name_;
}

IOdictionary& nestedBlendDsPimpleSolver::nestingControlDict()
{
    return nestingControlDict_;
}

Time& nestedBlendDsPimpleSolver::runTime()
{
    return runTime_;
}

fvMesh& nestedBlendDsPimpleSolver::mesh()
{
    return mesh_;
}

volVectorField& nestedBlendDsPimpleSolver::U()
{
    return U_;
}

volScalarField& nestedBlendDsPimpleSolver::p()
{
    return p_;
}

surfaceScalarField& nestedBlendDsPimpleSolver::phi()
{
    return phi_;
}

functionObjectList& nestedBlendDsPimpleSolver::functionObjects()
{
    return functionObjects_;
}

List< nestedBlendDsPimpleSolver* >& nestedBlendDsPimpleSolver::nestedSolvers()
{
    return nestedSolvers_;
}

bool nestedBlendDsPimpleSolver::runNested()
{
    return runNested_;
}

// * * * * * * * * * * * * * * * * set Functions * * * * * * * * * * * * * * //

void nestedBlendDsPimpleSolver::runNested(bool b)
{
    runNested_ = b;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
