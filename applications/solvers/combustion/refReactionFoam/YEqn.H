//if (!combustion->active())
{
    Info<<"# Species equations are being solved."<<endl;
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YEqn
            (
                fvm::ddt(rho, Yi)
              + scalarConvection->fvmDiv(phi, Yi)
                // Turbulent schmidt number is hard-coded at 1.0
              - fvm::laplacian(turbulence->muEff()/1.0, Yi)
//              - fvm::laplacian((
//                                 turbulence->mu()
//                                +turbulence->mut()/0.78
//                               ), Yi)
             ==
                parcels.SYi(i, Yi)
              //+ combustion->R(Yi)
              + fvOptions(rho, Yi)
            );

            YEqn.relax();

            fvOptions.constrain(YEqn);

            YEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
