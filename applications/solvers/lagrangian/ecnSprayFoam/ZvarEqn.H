{
    const scalar schmidt = 1.0;
    const scalar Cchi = 2.0;
    const volScalarField& k = turbulence->k();
    const volScalarField& e = turbulence->epsilon();
    const dimensionedScalar k_SMALL("k_SMALL", k.dimensions(), 1.e-09);

    Chi = Cchi*Zvar*e/(k+k_SMALL);

    fvScalarMatrix ZvarEqn
    (
         fvm::ddt(rho, Zvar)
       + scalarConvection->fvmDiv(phi, Zvar)
       - fvm::laplacian(turbulence->mut()/schmidt, Zvar)
//       - fvm::laplacian(turbulence->muEff()/combustion->schmidt(), Zvar)
//       - fvm::laplacian((
//                            turbulence->mu()
//                           +turbulence->mut()/combustion->schmidt()
//                        ), Zvar)
//       - 2.0*turbulence->mut()/combustion->schmidt()*magSqr(fvc::grad(fuelMFR))
      ==
         2.0*turbulence->mut()/schmidt*magSqr(fvc::grad(Z))
       - rho*Chi
       + fvOptions(rho, Zvar)
    );

    ZvarEqn.relax();

    fvOptions.constrain(ZvarEqn);

    ZvarEqn.solve(mesh.solver("Zvar"));

    fvOptions.correct(Zvar);

    Zvar.max(0.0);
}
