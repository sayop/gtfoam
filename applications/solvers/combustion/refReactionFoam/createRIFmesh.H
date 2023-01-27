Foam::Time rifRunTime
(
    Foam::Time::controlDictName,
    "RIFs",
    "rif"
);

//- set rifRunTime equivalent to CFD runTime. 
//- It matters when restart operates.
rifRunTime.setTime(runTime);
Info<<"timeName:"<<rifRunTime.timeName()<<endl;

Foam::fvMesh rifMesh
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        rifRunTime.timeName(),
        rifRunTime,
        Foam::IOobject::MUST_READ
    )
);

volScalarField I = rifMesh.C().component(0);
Info<<"sayop:"<<I<<endl;
Info<<"sayop:"<<rifRunTime.timeName()<<endl;
