label fuelIndex = -1;
forAll(Y, i)
{
    if (Y[i].name() == combustion->fuelSpecie())
    {
        fuelIndex = i;
    }
}

tmp<fv::convectionScheme<scalar>> scalarConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,scalar)")
    )
);


//- Update mixture fraction field as if it transports fuel specie with evaporation source.
fvScalarMatrix ZEqn
(
    fvm::ddt(rho, Z)
  + scalarConvection->fvmDiv(phi, Z)
  - fvm::laplacian(turbulence->muEff()/combustion->schmidt(), Z)
//  - fvm::laplacian((
//                       turbulence->mu()
//                      +turbulence->mut()/combustion->schmidt()
//                   ), Z)
 ==
    parcels.SYi(fuelIndex, Z)
  + fvOptions(rho, Z)
);

ZEqn.relax();

fvOptions.constrain(ZEqn);

ZEqn.solve(mesh.solver("Z"));

fvOptions.correct(Z);

Z.max(0.0);

volScalarField Ztot
(
    IOobject
    (
        "Ztot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Ztot", dimless, 0.0)
);

//- update the passive scalar equation for Zi: marker solution
forAll(Zi, i)
{
    volScalarField& Zii = Zi[i];
    //Info<<runTime.value()<<"&"<<combustion->tInitFlamelet(i)<<endl;
//    if (runTime.value() < combustion->tInitFlamelet(i))
//    {
//        continue;
//    }

    fvScalarMatrix ZiEqn
    (
        fvm::ddt(rho, Zii)
      + scalarConvection->fvmDiv(phi, Zii)
      - fvm::laplacian(turbulence->muEff()/combustion->schmidt(), Zii)
//      - fvm::laplacian((
//                           turbulence->mu()
//                          +turbulence->mut()/combustion->schmidt()
//                       ), Zii)
     ==
        parcels.SZYi(fuelIndex, i, Zii)
      + fvOptions(rho, Zii)
    );
 
    ZiEqn.relax();

    fvOptions.constrain(ZiEqn);

    ZiEqn.solve(mesh.solver("Zi"));

    fvOptions.correct(Zii);

    Zii.max(0.0);

    Ztot += Zii;
}


Ztot.max(1e-16);
// Adjust level of Zii to make Z bounded to 1 as maximum.
// This is to avoid errors that can be induced by numerically transporting 
// each of Zi quantity independently.
// When Z is transported, all Zi quantities are bounded by Z, but if Zi is 
// separately transported, applied numerical schemes (i.e. TVD method) are
// independently acting on each Zi equation, so any possible numerical errors
// is likely to make solution unbounded in Z solution.
forAll(Zi, i)
{
    volScalarField& Zii = Zi[i];
    Zii = Zii * Z / Ztot;
}
