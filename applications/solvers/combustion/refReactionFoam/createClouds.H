Info<< "\nConstructing reacting cloud" << endl;
//basicSprayCloud parcels
//(
//    "sprayCloud",
//    rho,
//    U,
//    g,
//    slgThermo
//);

basicSprayCloud parcels
(
    "sprayCloud",
    rho,
    U,
    Z,
    g,
    slgThermo,
    &combustion()
);

//const flameletModel* flamelet;
//flamelet = &combustion();
//Info<<"SAYOP: "<<flamelet->nFlamelets()<<endl;
//flamelet = combustion;
//flameletModel& flamelet = combustion();
//Foam::flameletModels::psiFlameletModel flamelet;
//autoPtr<flameletModels::psiFlameletModel> fRef(combustion);

if (parcels.solution().coupled() || parcels.injectors().size() >= 1)
{
    Info<<"Flamelet time history is set up based on first injector time history"<<endl;
    combustion->setFlameletTimes(parcels.injectors()[0].timeStart(),
                                 parcels.injectors()[0].timeEnd());
}

//Info<< "\nConstructing flamelet cloud" << endl;
//ameletSprayCloud flameletParcels
//
//  "flameletCloud",
//  rho,
//  U,
//  g,
//  slgThermo
//

//fo<< "SAYOP: Playing around with flamelet cloud" << endl;

//scalar gasMass0 = fvc::domainIntegrate(rho).value();
//Info<< "gasMass0=" << gasMass0 << endl;
