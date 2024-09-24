-- this is a dummy ideal gas gas model made by PITOT3 on the 24/06/2024 17:06:11
-- purely to perform a frozen shock over the test model.
-- Do not use it for anything else! As the entropy, viscosity and 
-- thermal conductivity values are just dummy values for air.

model = "IdealGas"

IdealGas = {
    speciesName = 'dummy PITOT3 test section species',
    mMass = 0.02885,
    gamma = 1.3934968622514174,
    entropyRefValues = {
        s1 = 0.0,
        T1 = 298.15,
        p1 = 101.325e3
    },
    viscosity = {
    model = 'Sutherland',
    mu_ref = 1.716e-5,
    T_ref = 273.0,
    S = 111.0,
    },
    thermCondModel = {
    model = 'Sutherland',
    T_ref = 273.0,
    k_ref = 0.0241,
    S = 194.0
    }
}

