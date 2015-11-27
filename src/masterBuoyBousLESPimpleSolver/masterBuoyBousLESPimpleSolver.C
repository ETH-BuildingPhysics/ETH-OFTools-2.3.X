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

#include "masterBuoyBousLESPimpleSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
//        defineTypeNameAndDebug(masterBuoyBousLESPimpleSolver, 0);
//        defineRunTimeSelectionTable(masterBuoyBousLESPimpleSolver, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

masterBuoyBousLESPimpleSolver::masterBuoyBousLESPimpleSolver
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
            Foam::fvMesh::defaultRegion,
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
    g_
    (
        IOobject
        (
            "g",
            runTime_.constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    T_
    (
        IOobject
        (
            "T",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    p_rgh_
    (
      IOobject
      (
          "p_rgh",
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
    beta_(laminarTransport_.lookup("beta")),
    TRef_(laminarTransport_.lookup("TRef")),
    Pr_(laminarTransport_.lookup("Pr")),
    Prt_(laminarTransport_.lookup("Prt")),
    turbulence_(incompressible::turbulenceModel::New(U_, phi_, laminarTransport_)),
    rhok_
    (
        IOobject
        (
            "rhok",
            runTime_.timeName(),
            mesh_
        ),
        1.0 - beta_*(T_ - TRef_)
    ),
    alphat_
    (
        IOobject
        (
            "alphat",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Rwall_
    (
        IOobject
        (
            "Rwall",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    qwall_
    (
        IOobject
        (
            "qwall",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    gh_("gh", g_ & mesh_.C()),
    ghf_("ghf", g_ & mesh_.Cf()),
    p_
    (
        IOobject
        (
            "p",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh_ + rhok_*gh_
    ),
    radiation_(radiation::radiationModel::New(T_)),
    rhoCpRef_("rhoCpRef",dimDensity*dimEnergy/dimMass/dimTemperature,1.0),
    fvOptions_(mesh_),
    pimple_(mesh_),
    pRefCell_(0),
    pRefValue_(0.0),
    cumulativeContErr_(0),
    nestedCases_(nestingControlDict_.lookup("nestedCases")),
    nestedSolvers_(),
    runNested_(runNested)
{
    Info << "Create mesh for time = " << runTime_.timeName() << nl << endl;
    Info << "Reading g\n" << endl;
    Info << "Reading thermophysical properties\n" << endl;
    Info << "Reading field T\n" << endl;
    Info << "Reading field p_rgh\n" << endl;
    Info << "Reading field U\n" << endl;
    Info << "Creating turbulence model\n" << endl;
    Info << "Reading field alphat\n" << endl;
    Info << "Reading field Rwall\n" << endl;
    Info << "Reading field qwall\n" << endl;
    Info << "Calculating field g.h\n" << endl;
    Info << "Reading/calculating face flux field phi\n" << endl;

    setRefCell
    (
        p_,
        p_rgh_,
        mesh_.solutionDict().subDict("PIMPLE"),
        pRefCell_,
        pRefValue_
    );
    if (p_rgh_.needReference())
    {
        p_ += dimensionedScalar
        (
            "p",
            p_.dimensions(),
            pRefValue_ - getRefCellValue(p_, pRefCell_)
        );
    };

    if (radiation_->radiation())
    {
        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                runTime_.constant(),
                runTime_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false  // do not register!
            )
        );
        dimensionedScalar rhoRef(transportProperties.lookup("rhoRef"));
        dimensionedScalar CpRef(transportProperties.lookup("CpRef"));
        rhoCpRef_ = rhoRef*CpRef;
    }

    functionObjects_.start();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

masterBuoyBousLESPimpleSolver::~masterBuoyBousLESPimpleSolver()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void masterBuoyBousLESPimpleSolver::setNestedSolvers
(
    List< nestedBlendDsBuoyBousLESPimpleSolver* > nestedSolvers
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


void masterBuoyBousLESPimpleSolver::runFunctionObjects()
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


void masterBuoyBousLESPimpleSolver::solve()
{
    Info << name_ << ": Time = " << runTime_.timeName() << nl << endl;

    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple_.loop())
    {
        #include "UEqn.H"
        #include "turbulenceCorrect.H"
        #include "TEqn.H"

        // --- Pressure corrector loop
        while (pimple_.correct())
        {
            #include "pEqn.H"
            #include "turbulenceCorrect.H"
            #include "TEqn.H"
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
}

// * * * * * * * * * * * * * * * Access Functions  * * * * * * * * * * * * * //
word& masterBuoyBousLESPimpleSolver::name()
{
    return name_;
}

Time& masterBuoyBousLESPimpleSolver::runTime()
{
    return runTime_;
}

IOdictionary& masterBuoyBousLESPimpleSolver::nestingControlDict()
{
    return nestingControlDict_;
}

fvMesh& masterBuoyBousLESPimpleSolver::mesh()
{
    return mesh_;
}

volScalarField& masterBuoyBousLESPimpleSolver::T()
{
    return T_;
}

volVectorField& masterBuoyBousLESPimpleSolver::U()
{
    return U_;
}

volScalarField& masterBuoyBousLESPimpleSolver::p_rgh()
{
    return p_rgh_;
}

volScalarField& masterBuoyBousLESPimpleSolver::p()
{
    return p_;
}

surfaceScalarField& masterBuoyBousLESPimpleSolver::phi()
{
    return phi_;
}

functionObjectList& masterBuoyBousLESPimpleSolver::functionObjects()
{
    return functionObjects_;
}

List< nestedBlendDsBuoyBousLESPimpleSolver* >& masterBuoyBousLESPimpleSolver::nestedSolvers()
{
    return nestedSolvers_;
}

bool masterBuoyBousLESPimpleSolver::runNested()
{
    return runNested_;
}

// * * * * * * * * * * * * * * * * set Functions * * * * * * * * * * * * * * //

void masterBuoyBousLESPimpleSolver::runNested(bool b)
{
    runNested_ = b;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
