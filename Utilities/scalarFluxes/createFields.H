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



    wordList names(transportProperties.lookup("scalars"));

    PtrList<volScalarField> sMean(names.size());


    
    Info<< "Creating Mean Scalars\n" << endl;
    forAll(sMean, i)
    {
        IOobject header
        (
            names[i]+EXT_MEAN,
            "MeanFlow",
            mesh,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.headerOk())
        {
            sMean.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        names[i]+EXT_MEAN,
                        "MeanFlow",
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                )
            );
        }
    }



    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );
