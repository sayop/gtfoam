{
    volScalarField& he = thermo.he();
//    solve
//    (
//        fvm::ddt(rho, he) - chemistry.Sh(),
//        mesh.solver("h")
//    );
    fvScalarMatrix heEqn
    (
        rho*invRhoChiSt*fvm::ddt(he)
//        rho*fvm::ddt(he)
//      - fvm::laplacian(0.5*rhoChiSt, he)
      - fvm::laplacian(he)
     ==
        invRhoChiSt*chemistry.Sh()
    );
    heEqn.relax();
    heEqn.solve();

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
