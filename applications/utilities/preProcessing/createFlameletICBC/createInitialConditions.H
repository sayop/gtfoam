    Info<< "Reading initial conditions.\n" << endl;
    IOdictionary initialConditions
    (
        IOobject
        (
            "initialConditions",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalar p0 = readScalar(initialConditions.lookup("p"));
    scalar Tair = readScalar(initialConditions.lookup("Tair"));
    scalar Tfuel = readScalar(initialConditions.lookup("Tfuel"));

    #include "createTemporaryFields.H"
    #include "createThermoFields.H"

    word constProp(initialConditions.lookup("constantProperty"));
    if ((constProp != "pressure") && (constProp != "volume"))
    {
        FatalError << "in initialConditions, unknown constantProperty type "
            << constProp << nl << " Valid types are: pressure volume."
            << abort(FatalError);
    }

    word fractionBasis(initialConditions.lookup("fractionBasis"));
    if ((fractionBasis != "mass") && (fractionBasis != "mole"))
    {
        FatalError << "in initialConditions, unknown fractionBasis type " << nl
            << "Valid types are: mass or mole."
            << fractionBasis << abort(FatalError);
    }


    label nSpecie = Y.size();
    PtrList<gasHThermoPhysics> specieData(nSpecie);
    forAll(specieData, i)
    {
        specieData.set
        (
            i,
            new gasHThermoPhysics
            (
                dynamic_cast<const reactingMixture<gasHThermoPhysics>&>
                    (thermo).speciesData()[i]
            )
        );
    }



    scalarList Y0air(nSpecie, 0.0);
    scalarList Y0fuel(nSpecie, 0.0);
    scalarList X0air(nSpecie, 0.0);    
    scalarList X0fuel(nSpecie, 0.0);    

    dictionary airFractions(initialConditions.subDict("air"));
    dictionary fuelFractions(initialConditions.subDict("fuel"));

    if (fractionBasis == "mole")
    {
        forAll(Y, i)
        {
            const word& name = Y[i].name();
            if (airFractions.found(name))
            {
                X0air[i] = readScalar(airFractions.lookup(name));
            }
            if (fuelFractions.found(name))
            {
                X0fuel[i] = readScalar(fuelFractions.lookup(name));
            }
        }

        scalar mwAir = 0.0;
        scalar mwFuel = 0.0;
        const scalar mTotAir = sum(X0air);
        const scalar mTotFuel = sum(X0fuel);
        forAll(Y, i)
        {
            X0air[i] /= mTotAir;
            X0fuel[i] /= mTotFuel;
            mwAir += specieData[i].W()*X0air[i];
            mwFuel += specieData[i].W()*X0fuel[i];
        }

        forAll(Y, i)
        {
            Y0air[i] = X0air[i]*specieData[i].W()/mwAir;
            Y0fuel[i] = X0fuel[i]*specieData[i].W()/mwFuel;
        }
    }
    else  // mass fraction
    {
        forAll(Y, i)
        {
            const word& name = Y[i].name();
            if (airFractions.found(name))
            {
                Y0air[i] = readScalar(airFractions.lookup(name));
            }
            if (fuelFractions.found(name))
            {
                Y0fuel[i] = readScalar(fuelFractions.lookup(name));
            }
        }

        scalar invWair = 0.0;
        scalar invWfuel = 0.0;
        const scalar mTotAir = sum(Y0air);
        const scalar mTotFuel = sum(Y0fuel);
        forAll(Y, i)
        {
            Y0air[i] /= mTotAir;
            Y0fuel[i] /= mTotFuel;
            invWair += Y0air[i]/specieData[i].W();
            invWfuel += Y0fuel[i]/specieData[i].W();
        }
        const scalar mwAir = 1.0/invWair;
        const scalar mwFuel = 1.0/invWfuel;

        forAll(Y, i)
        {
            X0air[i] = Y0air[i]*mwAir/specieData[i].W();
            X0fuel[i] = Y0fuel[i]*mwFuel/specieData[i].W();
        }
    }

    forAll(Y, i)
    {
        Y[i] = (1 - Z) * Y0air[i] + Z * Y0fuel[i];
    }

    volScalarField Tboundary = (1 - Z) * Tair + Z * Tfuel;
    Tboundary.dimensions().reset(dimTemperature);
    T = Tboundary;

    scalar hAir = 0.0;
    scalar hFuel = 0.0;

    forAll(Y, i)
    {
        hAir += Y0air[i] * specieData[i].Hs(p0, Tair);
        hFuel += Y0fuel[i] * specieData[i].Hs(p0, Tfuel);
    }
    //- Since he() is dimensioned scalar, =operation can't be used for he() along with Z.
    volScalarField hMix = (1 - Z) * hAir + Z * hFuel;
    hMix.dimensions().reset(dimEnergy/dimMass);
    thermo.he() = hMix;

//    thermo.he().write();
    thermo.correct();

    //- Re-write field data for correct initial condition
    forAll(Y, i)
    {
        volScalarField Yi
        (
            IOobject
            (
                Y[i].name(),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Y[i],
            "fixedValue"
        );
        Yi.write();
    }

    volScalarField TfixedBC
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T,
        "fixedValue"
    );
    TfixedBC.write();
