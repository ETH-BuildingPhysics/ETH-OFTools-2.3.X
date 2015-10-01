    Info<< "\nReading field UMean\n" << endl;
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
    
    wordList names;
    
    if (args.optionFound("scalars")==true)
    {
        names = args.optionReadList<word>("scalars");
    }
    else
    {
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

        wordList tmp(transportProperties.lookup("scalars"));
        names = tmp;
    }

    PtrList<volScalarField> sMean(names.size());
    PtrList<volVectorField> usMean(names.size());

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
            Info<< "Reading field " << names[i]+EXT_MEAN << endl;
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

    Info << endl;

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
            Info<< "Reading field " << "U"+names[i]+EXT_MEAN << endl;
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

    Info << endl;
