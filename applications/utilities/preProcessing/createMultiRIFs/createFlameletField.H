    Info<< "Reading flamelet setups.\n" << endl;
    IOdictionary flameletConditions
    (
        IOobject
        (
            "flameletProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );


    label nFlamelets = readLabel(flameletConditions.lookup("nFlamelets"));

    const fileName solutionDir = cwd() + "/RIFs";
    const fileName caseRepoDir = cwd()+"/RIFs/rif";
    List<fileName> caseDir;
    caseDir.setSize(nFlamelets);

    IFstream is(solutionDir);
    if (!is.good())
    {
        mkDir(solutionDir);
        FatalErrorInFunction
            << "Failed opening flamelet solution directory! Just created directory: "
            << is.name() << endl
            << "Please make a reference RIF case setup in the directory!!"
            << exit(FatalError);
    }

    IFstream isCaseRepo(caseRepoDir);
    if (!isCaseRepo.good())
    {
        mkDir(caseRepoDir);
        FatalErrorInFunction
            << "Failed reading flamelet case file, "
            << isCaseRepo.name() << endl
            << "Please make a reference case setup for the flamelet solver!!"
            << exit(FatalError);
    }

    if (forceOverwrite)
    {
        Info<< "Removing all case directory..." << endl;
        system("rm -rf "+solutionDir+"/rif.*");
    }

    forAll (caseDir, i)
    {
        caseDir[i] = solutionDir/"rif."+std::to_string(i);
        IFstream is(caseDir[i]);
        if (!is.good())
        {
            Info<<"Creating directory: "<<caseDir[i]<<endl;
            mkDir(caseDir[i]);
            system("cp -r "+caseRepoDir+"/0 "+caseDir[i]+"/");
            system("cp -r "+caseRepoDir+"/constant "+caseDir[i]+"/");
            system("cp -r "+caseRepoDir+"/system "+caseDir[i]+"/");
            system("cp -r "+caseRepoDir+"/chemkin "+caseDir[i]+"/");
            IFstream isPar(caseRepoDir+"/processor0");
            if (isPar.good())
            {
                system("cp -r "+caseRepoDir+"/processor* "+caseDir[i]+"/");
            }
        }
        else
        {
            Info<<caseDir[i]<<" already exists!"<<endl;
        }
  
    }



