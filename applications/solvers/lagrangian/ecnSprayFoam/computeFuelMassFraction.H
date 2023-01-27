{
    scalarField liquidMass(Z.size(), 0.0);
    scalarField vaporMass = Z * rho * mesh.V();
    scalarField totalMass = rho * mesh.V();

//    label fuelIndex = -1;
//    forAll(Y, i)
//    {
//        if (Y[i].name() == fuelSpecie)
//        {
//            fuelIndex = i;
//        }
//    }

//    volScalarField& vap  = Y[fuelIndex];
//    scalarField liquidMass(vap.size(), 0.0);
//    scalarField vaporMass = vap * rho * mesh.V();
//    scalarField totalMass = rho * mesh.V();

    tmp<volScalarField> tMFR
    (
        new volScalarField
        (
            IOobject
            (
                "tMFR",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    forAllIter(basicSprayCloud, parcels, pIter)
    {
        liquidMass[pIter().cell()] += pIter().mass() * pIter().nParticle();
    }
    totalMass += liquidMass;

    forAll(tMFR.ref(), celli)
    {
        tMFR.ref()[celli] = (liquidMass[celli]+vaporMass[celli])/totalMass[celli];
    }

    fuelMFR = tMFR;
    fuelMFR.min(1.0);
}
