psi = thermo.psi();

Info<< "min/max(psi) = "
    << min(psi).value() << ", " << max(psi).value() << endl;

p = totalMass0/fvc::domainIntegrate(psi); // p0 update
p.correctBoundaryConditions();

Info<< "min/max(p) = "
    << min(p).value() << ", " << max(p).value() << endl;

rho = thermo.rho(); // function of T only so no need to wait until p0 update

dimensionedScalar totalMass = fvc::domainIntegrate(rho);
scalar constantMass = ((scalar(1.0)- totalMass/totalMass0)*100.0).value();

Info<< "Percentage mass loss/gain  " << constantMass << endl;

    Info<< "min/max(rho) = "
        << min(rho).value() << ", " << max(rho).value() << endl;

    Info<< "min/max(cellLength) = "
        << min(cellLength).value() << ", " << max(cellLength).value() << endl;

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA("HbyA", U); 
HbyA = rAU*UEqn.H(); 
surfaceScalarField phig("phig", -rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(rho*HbyA)
      + rhorAUf*fvc::ddtCorr(rho, U, phi)
    )
    + phig
);

constrainPressure(pd, rho, U, phiHbyA, rhorAUf);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pdEqn
    (
        fvm::ddt(psi, pd)
      - fvc::ddt(psi, pd)
      + fvc::ddt(rho) // what if this should be ddt(psi)*p
      + fvc::div(phiHbyA)
      - fvm::laplacian(rhorAUf, pd)
     ==
        fvOptions(psi, pd, rho.name())
    );

    pdEqn.solve(mesh.solver(pd.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA + pdEqn.flux();
        U = HbyA + rAU*fvc::reconstruct((pdEqn.flux() + phig)/rhorAUf);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}


