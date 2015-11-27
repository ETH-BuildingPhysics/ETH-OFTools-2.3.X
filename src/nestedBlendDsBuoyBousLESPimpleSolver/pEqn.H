{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U_);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phig(-rAUf*ghf_*fvc::snGrad(rhok_)*mesh_.magSf());

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh_.Sf())
      + rAUf*fvc::ddtCorr(U_, phi_)
      + phig
    );

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh_.boundaryField(),
        (
            phiHbyA.boundaryField()
          - (mesh_.Sf().boundaryField() & U_.boundaryField())
        )/(mesh_.magSf().boundaryField()*rAUf.boundaryField())
    );

    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh_) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell_, getRefCellValue(p_rgh_, pRefCell_));

        p_rghEqn.solve(mesh_.solver(p_rgh_.select(pimple_.finalInnerIter())));

        if (pimple_.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi_ = phiHbyA - p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh_.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U_ = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U_.correctBoundaryConditions();
            fvOptions_.correct(U_);
        }
    }

    #include "continuityErrs.H"

    p_ = p_rgh_ + rhok_*gh_;

    if (p_rgh_.needReference())
    {
        p_ += dimensionedScalar
        (
            "p",
            p_.dimensions(),
            pRefValue_ - getRefCellValue(p_, pRefCell_)
        );
        p_rgh_ = p_ - rhok_*gh_;
    }
}
