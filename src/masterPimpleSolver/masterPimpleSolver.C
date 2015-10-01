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

#include "masterPimpleSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
//        defineTypeNameAndDebug(masterPimpleSolver, 0);
//        defineRunTimeSelectionTable(masterPimpleSolver, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

masterPimpleSolver::masterPimpleSolver
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
    runNested_(runNested)
{
    Info << "Create mesh for time = " << runTime_.timeName() << nl << endl;
    Info << "Reading field p\n" << endl;
    Info << "Reading field U\n" << endl;
    Info << "Reading/calculating face flux field phi\n" << endl;

    setRefCell(p_, mesh_.solutionDict().subDict("PIMPLE"), pRefCell_, pRefValue_);
    functionObjects_.start();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

masterPimpleSolver::~masterPimpleSolver()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void masterPimpleSolver::setNestedSolvers
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


void masterPimpleSolver::runFunctionObjects()
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


void masterPimpleSolver::solve()
{
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
}

// * * * * * * * * * * * * * * * Access Functions  * * * * * * * * * * * * * //
word& masterPimpleSolver::name()
{
    return name_;
}

Time& masterPimpleSolver::runTime()
{
    return runTime_;
}

IOdictionary& masterPimpleSolver::nestingControlDict()
{
    return nestingControlDict_;
}

fvMesh& masterPimpleSolver::mesh()
{
    return mesh_;
}

volVectorField& masterPimpleSolver::U()
{
    return U_;
}

volScalarField& masterPimpleSolver::p()
{
    return p_;
}

surfaceScalarField& masterPimpleSolver::phi()
{
    return phi_;
}

functionObjectList& masterPimpleSolver::functionObjects()
{
    return functionObjects_;
}

List< nestedBlendDsPimpleSolver* >& masterPimpleSolver::nestedSolvers()
{
    return nestedSolvers_;
}

bool masterPimpleSolver::runNested()
{
    return runNested_;
}

// * * * * * * * * * * * * * * * * set Functions * * * * * * * * * * * * * * //

void masterPimpleSolver::runNested(bool b)
{
    runNested_ = b;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
