// Solve the Momentum equation

fvVectorMatrix UEqn
(
    fvm::ddt(U_)
  + fvm::div(phi_, U_)
  + turbulence_->divDevReff(U_)
  + fvc::div(Rwall_)
 ==
    fvc::Sp(wDs_/DsRelaxTime,UDsBlended_)
  - fvc::Sp(wDs_/DsRelaxTime,U_)
  + fvOptions_(U_)
);

UEqn.relax();

fvOptions_.constrain(UEqn);

if (pimple_.momentumPredictor())
{
    Foam::solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
              - ghf_*fvc::snGrad(rhok_)
              - fvc::snGrad(p_rgh_)
            )*mesh_.magSf()
        )
    );

    fvOptions_.correct(U_);
}
