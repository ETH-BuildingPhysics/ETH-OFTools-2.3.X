surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

volVectorField HbyA("HbyA", U_);
HbyA = rAU*UEqn().H();

if (pimple_.nCorrPISO() <= 1)
{
    UEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (fvc::interpolate(HbyA) & mesh_.Sf())
  + rAUf*fvc::ddtCorr(U_, phi_)
);

fvOptions_.makeRelative(phiHbyA);

adjustPhi(phiHbyA, U_, p_);

// Non-orthogonal pressure corrector loop
while (pimple_.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUf, p_) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell_, pRefValue_);

    pEqn.solve(mesh_.solver(p_.select(pimple_.finalInnerIter())));

    if (pimple_.finalNonOrthogonalIter())
    {
        phi_ = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p_.relax();

U_ = HbyA - rAU*fvc::grad(p_);
U_.correctBoundaryConditions();
fvOptions_.correct(U_);
