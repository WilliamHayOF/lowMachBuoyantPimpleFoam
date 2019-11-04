Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));
rhoThermo& thermo = pThermo();
thermo.validate(args.executable(), "h");

volScalarField& T = thermo.T();
T.correctBoundaryConditions();

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
rho.oldTime(); // need to store rho for ddt term ? like in multiphase ? 

dimensionedScalar totalMass0 = fvc::domainIntegrate(rho); // total mass based on thermo.rho AT t0

volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh
    ),
    thermo.psi()
);
psi.oldTime(); // need to store for ddt term regardless of if we use ddt(rho) or ddt(psi)*p

volScalarField& p = thermo.p();

p = totalMass0/fvc::domainIntegrate(psi);
p.correctBoundaryConditions();

Info<< "Reading field U\n" << endl;
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

#include "compressibleCreatePhi.H"

Info<< "Creating field pd\n" << endl;
volScalarField pd
(
    IOobject
    (
        "pd",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt //dummy dpdt kept to keep ternary operator happy (EEqn.H)
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

mesh.setFluxRequired(pd.name());

#include "createFvOptions.H"


