    // Thermal expansion coefficient [1/K]
    dimensionedScalar Pr(laminarTransport.lookup("Pr"));

    // Turbulent Prandtl number
    dimensionedScalar Prt(laminarTransport.lookup("Prt"));

    // Heat capacity
    dimensionedScalar Cp(laminarTransport.lookup("Cp"));

    // Fluid density
    dimensionedScalar rho(laminarTransport.lookup("rho"));
