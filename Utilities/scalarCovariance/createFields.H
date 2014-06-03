    Info<< "Reading field UMean\n" << endl;
    volVectorField UMean
    (
        IOobject
        (
            "UMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

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
    PtrList<volVectorField> usMean(names.size());

    
    Info<< "Creating Mean Scalars\n" << endl;
    forAll(sMean, i)
    {
        IOobject header
        (
            names[i]+EXT_MEAN,
            runTime.timeName(),
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
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                )
            );
        }
    }

    forAll(usMean, i)
    {
        IOobject header
        (
            "U"+names[i]+EXT_MEAN,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.headerOk())
        {
            usMean.set
            (
                i,
                new volVectorField
                (
                    IOobject
                    (
                        "U"+names[i]+EXT_MEAN,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                )
            );
        }
    }
