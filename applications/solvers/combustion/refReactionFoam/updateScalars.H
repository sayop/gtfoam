//- update CFD species field from RIF solution.
if (combustion->updateSpecies(Y))
{
    Info<<":-) Species field updated!"<<endl;
}

//- update CFD energy field from RIF solution
volScalarField& he = thermo.he();
if (combustion->updateEnthalpy(he))
{
    Info<<":-) Enthalpy field updated!"<<endl;
}

thermo.correct();
Info<<nl<< "# T gas min/max   = " << min(T).value() << ", "
    << max(T).value() << endl;

//- compute heat release quantity
dQ = combustion->dQ();
