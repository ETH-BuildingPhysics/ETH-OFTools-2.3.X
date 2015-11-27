/*---------------------------------------------------------------------------*\
                         create vol***Field
                      (can be plotted in paraview)
\*---------------------------------------------------------------------------*/

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

    
Info << "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info << "Creating field alphaEff\n" << endl;
volScalarField alphaEff
(
    IOobject
    (
        "alphaEff",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alphaEff", dimLength*dimLength/dimTime,0)
);


Info << "Creating field alphat\n" << endl;
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alphat", dimLength*dimLength/dimTime,0)
);


/*---------------------------------------------------------------------------*/

    
#   include "createPhi.H"

//    label pRefCell = 0;
//    scalar pRefValue = 0.0;
//    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    turbulence->read();
