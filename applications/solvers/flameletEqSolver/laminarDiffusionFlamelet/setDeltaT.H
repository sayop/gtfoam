if (adjustTimeStep)
{
    scalarList dtChemList(Pstream::nProcs(),0.0);
    dtChemList[Pstream::myProcNo()] = dtChem;
    // Gather data onto master
    Pstream::gatherList(dtChemList);

    dtChem = min(dtChemList);
    Pstream::scatter(dtChem);

    runTime.setDeltaT(min(dtChem, maxDeltaT));
    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
}

