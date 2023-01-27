label fuelIndex = -1;
forAll(Y, i)
{
    if (Y[i].name() == fuelSpecie)
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
  - fvm::laplacian(turbulence->muEff()/1.0, Z)
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
