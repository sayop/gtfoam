    Info<< nl << "Reading thermophysicalProperties" << endl;
    autoPtr<psiChemistryModel> pChemistry(psiChemistryModel::New(mesh));

    psiChemistryModel& chemistry = pChemistry();

    psiReactionThermo& thermo = chemistry.thermo();
    basicMultiComponentMixture& composition = thermo.composition();

    thermo.validate(args.executable(), "h");

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    Info<< "Reading mixing conditions.\n" << endl;
    IOdictionary mixingConditions
    (
        IOobject
        (
            "mixingConditions",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalar scalarChiSt = max(SMALL,readScalar(mixingConditions.lookup("chiSt")));

    volScalarField chiSt
    (
        IOobject
        (
            "chiSt",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("chiSt", dimless/dimTime, scalarChiSt)
    );


//    volVectorField U
//    (
//        IOobject
//        (
//            "U",
//            runTime.timeName(),
//            runTime,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh,
//        dimensionedVector("zero", dimVelocity, Zero)
//    );

    volScalarField Z = mesh.C().component(0);
    Z.dimensions().reset(dimless);

//    #include "createPhi.H"

    PtrList<volScalarField>& Y = composition.Y();
//    volScalarField& p = thermo.p();
    volScalarField& T = thermo.T();

    scalar dtChem = refCast<const psiChemistryModel>(chemistry).deltaTChem()[0];


//
//- Following lines are defined for data collection across parallel runs.
//
    //- nDataSize: number of bins in mixture fraction
    label nDataSize = Z.size();

    //- When parallel running, nDataSize on each of processor IS NOT representative of
    //- total size of data. Thus, The actual data size should be sum of nDataSize from
    //- all processors.
    if (Pstream::parRun())
    {
        scalarList nDataSizeList(Pstream::nProcs(),0.0);
        nDataSizeList[Pstream::myProcNo()] = nDataSize;
        Pstream::gatherList(nDataSizeList);
        if (Pstream::master())
        {
            nDataSize = sum(nDataSizeList);
        }
        Pstream::scatter(nDataSize);
    }

    label nSpecie = Y.size();
    List<List<scalar>> massFraction(nSpecie);
    forAll(massFraction, speci)
    {
        massFraction[speci].setSize(nDataSize,0.0);
    }
    List<scalar> temperature(nDataSize,0.0);

    SortableList<scalar> mixtureFrac(nDataSize,0.0);
    List<label> startIndexList(Pstream::nProcs(),0);
    List<label> sortedIndx(nDataSize,0.0);

    //- Gather decomposed mixture fraction domain field.
    if (Pstream::parRun())
    {
        List<List<scalar>> procMixtureFrac(Pstream::nProcs());
        List<scalar>& myMixtureFrac = procMixtureFrac[Pstream::myProcNo()];
        myMixtureFrac.setSize(Z.size(), 0.0);
        forAll(Z, celli)
        {
            myMixtureFrac[celli] = Z[celli];
        }
        Pstream::gatherList(procMixtureFrac);

        if (Pstream::master())
        {
            label startIndex = 0;
            for (label proci = 0; proci < Pstream::nProcs(); proci++)
            {
                label dataSize = procMixtureFrac[proci].size();
                startIndexList[proci] = startIndex;
                SubList<scalar>
                (
                    mixtureFrac,
                    dataSize,
                    startIndex
                ) = procMixtureFrac[proci];
                startIndex += dataSize;
            }  
        }
        Pstream::scatter(startIndexList);
        // sort mixtureFrac list into ascending order
        mixtureFrac.sort();
        sortedIndx = mixtureFrac.indices();
        Pstream::scatter(mixtureFrac);        
        Pstream::scatter(sortedIndx);        
    }
    else
    {
        forAll(Z, celli)
        {
            mixtureFrac[celli] = Z[celli];
        }
    }

    //- Print out 1D flamelet data field in separate files
    PtrList<OFstream> specieDataFilePtrs_;
    specieDataFilePtrs_.setSize(nSpecie);
    fileName outputFileDir(args.rootPath()/args.globalCaseName()/"flameletDataOutput");

    mkDir(outputFileDir);
    OFstream fileTemperature(outputFileDir/"temperature_"+runTime.timeName()+".csv");

    if (Pstream::master())
    {
        fileTemperature<<runTime.timeName()<<",";
        forAll(mixtureFrac, celli)
        {
            fileTemperature<<mixtureFrac[celli]<<",";
        }
        fileTemperature<<endl;

        forAll(Y, speci)
        {
            word fName(Y[speci].name()+"_"+runTime.timeName()+".csv");
            specieDataFilePtrs_.set(speci, new OFstream(outputFileDir/fName));
            OFstream& file = specieDataFilePtrs_[speci];
            file<<runTime.timeName()<<",";    
            forAll(mixtureFrac, celli)
            {
                file<<mixtureFrac[celli]<<",";
            }
            file<<endl;
        }
    }
