{
    alphat_ = turbulence_->nut()/Prt_;
    alphat_.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence_->nu()/Pr_ + alphat_);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T_)
      + fvm::div(phi_, T_)
      - fvm::laplacian(alphaEff, T_)
      - fvc::div(qwall_)
     ==
        radiation_->ST(rhoCpRef_, T_)
      + fvOptions_(T_)
    );

    TEqn.relax();

    fvOptions_.constrain(TEqn);

    TEqn.solve();

    radiation_->correct();

    fvOptions_.correct(T_);

    rhok_ = 1.0 - beta_*(T_ - TRef_);
}
