tmp<Foam::volScalarField> tFnZ
(
    new volScalarField
    (
        IOobject
        (
            "fnZ",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Z
    )
);
volScalarField& fnZ = tFnZ.ref();
volScalarField::Boundary& Bf = fnZ.boundaryFieldRef();

forAll(Bf, patchi)
{
    forAll(Bf[patchi], facei)
    {
        scalar& value = Bf[patchi][facei];
        value = std::exp(-2.0*std::pow(erfc_inv(min(max(2.0*value,SMALL),1.9999999)),2.0));
    }
}
forAll(fnZ, celli)
{
    scalar& value = fnZ[celli];
    value = std::exp(-2.0*std::pow(erfc_inv(min(max(2.0*value,SMALL),1.9999999)),2.0));
}
scalar fZst = std::exp(-2.0*std::pow(erfc_inv(min(max(2.0*0.04509,SMALL),1.9999999)),2.0));

volScalarField ChiZ = fnZ * scalarChiSt / fZst;

volScalarField rhoChiSt
(
    IOobject
    (
        "rhoChiSt",
        runTime.timeName(),
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    0.5*chiSt*rho
);

volScalarField invRhoChiSt
(
    IOobject
    (
        "invRhoChiSt",
        runTime.timeName(),
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    1.0/rhoChiSt
);

//dimensionedScalar unity = dimensionedScalar("unity", dimless, 1.0);


//- Reset dimension for the purpose of laplacian term calculation 
//- in transport equation solver because in actual calculation, 
//- Z is considered length dimensioned. Therefore, in laplacian term, 
//- the length dimension should be canceled out.
rhoChiSt.dimensions().reset(dimMass/dimLength/dimTime);
invRhoChiSt.dimensions().reset(dimLength*dimTime/dimMass);

{
    forAll(Y, specieI)
    {
        volScalarField& Yi = Y[specieI];

//        solve
//        (
//            fvm::ddt(rho, Yi) - chemistry.RR(specieI),
//            mesh.solver("Yi")
//        );
        fvScalarMatrix YiEqn
        (
            rho*invRhoChiSt*fvm::ddt(Yi)
//          - fvm::laplacian(0.5*rhoChiSt, Yi)
//          + fvc::laplacian(0.5*rhoChiSt, Yi)
          - fvm::laplacian(Yi)
         ==
            invRhoChiSt*chemistry.RR(specieI)
        );
        YiEqn.relax();
        YiEqn.solve(mesh.solver("Yi"));
    }
}
