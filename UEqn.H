    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
      fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn 
         == 
            fvc::reconstruct
	    (
               (
		   - ghf*fvc::snGrad(rho)
                   -fvc::snGrad(pd)
	       )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
