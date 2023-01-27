#include "readGravitationalAcceleration.H"

Info<< "Creating flamelet model\n" << endl;

autoPtr<flameletModels::psiFlameletModel> combustion
(
    flameletModels::psiFlameletModel::New(mesh)
);

psiReactionThermo& thermo = combustion->thermo();
thermo.validate(args.executable(), "h", "e");


basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));

if (!composition.contains(inertSpecie))
{
    FatalErrorInFunction
        << "Specified inert specie '" << inertSpecie << "' not found in "
        << "species list. Available species:" << composition.species()
        << exit(FatalError);
}

volScalarField& p = thermo.p();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
