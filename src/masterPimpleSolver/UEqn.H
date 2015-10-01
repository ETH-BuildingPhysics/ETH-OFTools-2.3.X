// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U_)
  + fvm::div(phi_, U_)
  + turbulence_->divDevReff(U_)
 ==
    fvOptions_(U_)
);

UEqn().relax();

fvOptions_.constrain(UEqn());

volScalarField rAU(1.0/UEqn().A());

if (pimple_.momentumPredictor())
{
    Foam::solve(UEqn() == -fvc::grad(p_));

    fvOptions_.correct(U_);
}
