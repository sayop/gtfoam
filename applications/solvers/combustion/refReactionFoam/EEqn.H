tmp<fv::convectionScheme<scalar>> energyConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,energy)")
    )
);

{
    Info<<"# Energy equation is being solved."<<endl;
    volScalarField& he = thermo.he();

    // enthalpy balance equation
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + energyConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        rho*(U&g)
      + parcels.Sh(he)
      + radiation->Sh(thermo)
      + combustion->Sh()
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    radiation->correct();

    //Info<< "T gas min/max   = " << min(T).value() << ", "
    //    << max(T).value() << endl;
}
