fvScalarMatrix ZvarEqn
(
    fvm::ddt(rho, Zvar)
  + scalarConvection->fvmDiv(phi, Zvar)
  - fvm::laplacian(turbulence->mut()/combustion->schmidt(), Zvar)
//  - fvm::laplacian(turbulence->muEff()/combustion->schmidt(), Zvar)
//  - fvm::laplacian((
//                       turbulence->mu()
//                      +turbulence->mut()/combustion->schmidt()
//                   ), Zvar)
//  - 2.0*turbulence->mut()/combustion->schmidt()*magSqr(fvc::grad(fuelMFR))
 ==
    2.0*turbulence->mut()/combustion->schmidt()*magSqr(fvc::grad(Z))
  - rho*Chi
  + fvOptions(rho, Zvar)
);

ZvarEqn.relax();

fvOptions.constrain(ZvarEqn);

ZvarEqn.solve(mesh.solver("Zvar"));

fvOptions.correct(Zvar); 

Zvar.max(0.0);
