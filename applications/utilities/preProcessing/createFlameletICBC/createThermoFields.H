    Info<< nl << "Reading thermophysicalProperties" << endl;
    autoPtr<psiChemistryModel> pChemistry(psiChemistryModel::New(mesh));

    psiChemistryModel& chemistry = pChemistry();

    psiReactionThermo& thermo = chemistry.thermo();
    basicMultiComponentMixture& composition = thermo.composition();

    PtrList<volScalarField>& Y = composition.Y();
    volScalarField& T = thermo.T();


    volScalarField Z = mesh.C().component(0);
    Z.dimensions().reset(dimless);


