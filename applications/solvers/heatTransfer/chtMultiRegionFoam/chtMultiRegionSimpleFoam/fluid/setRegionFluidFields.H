    const fvMesh& mesh = fluidRegions[i];

    basicRhoThermo& thermo = thermoFluid[i];
    volScalarField& rho = rhoFluid[i];
    volScalarField& K = KFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    compressible::turbulenceModel& turb = turbulence[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();
    volScalarField& h = thermo.h();

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    radiation::radiationModel& rad = radiation[i];

    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];

    volScalarField& p_rgh = p_rghFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];
