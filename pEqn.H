psi = thermo.psi();

Info<< "min/max(psi) = "
    << min(psi).value() << ", " << max(psi).value() << endl;

//p = totalMass0/fvc::domainIntegrate(psi);
//p.correctBoundaryConditions();

Info<< "min/max(p) = "
    << min(p).value() << ", " << max(p).value() << endl;

rho = thermo.rho();

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
      + fvc::ddt(rho)
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


