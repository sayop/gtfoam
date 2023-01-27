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

fvScalarMatrix ZEqn
(
    fvm::ddt(rho, Z)
  + scalarConvection->fvmDiv(phi, Z)
  - fvm::laplacian(turbulence->muEff()/combustion->schmidt(), Z)
 ==
    parcels.SYi(fuelIndex, Z)
  + fvOptions(rho, Z)
);

ZEqn.relax();

fvOptions.constrain(ZEqn);

ZEqn.solve(mesh.solver("Z"));

fvOptions.correct(Z); 

Z.max(0.0);
Z.min(1.0);
