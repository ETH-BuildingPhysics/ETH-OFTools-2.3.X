    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
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

#   include "createPhi.H"

    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar Ds
    (
         transportProperties.lookup("Ds")
    );

    // laminar prandtl number [Pr]
    dimensionedScalar Pr(transportProperties.lookup("Pr"));
    //turbulent Prandtl number [Prt]
    dimensionedScalar Prt(transportProperties.lookup("Prt"));

    Info<< "Creating Scalars\n" << endl;


    wordList names(transportProperties.lookup("scalars"));
    PtrList<volScalarField> s(names.size());

    forAll(s, i)
    {
        IOobject header
        (
            names[i],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.headerOk())
        {
            s.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        names[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            volScalarField sdefault
            (
                IOobject
                (
                    "sdefault",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            s.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        names[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    sdefault
                )
            );
        }
    }
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    forAll(s, i)
	{
	    fields.add(s[i]);
	}

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    IObasicSourceList sources(mesh);
