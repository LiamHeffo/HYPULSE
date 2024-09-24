-- modified gas model made by PITOT3 on the 21/06/2024 17:59:49 without ions for low temperature operation

model = "CEAGas"

CEAGas = {
  mixtureName = "air13species-without-ions",
  speciesList = {'N2', 'O2', 'Ar', 'N', 'O', 'NO'},
  reactants = {N2 = 0.7811, O2 = 0.2095, Ar = 0.0093},
  inputUnits = 'moles',
  withIons = false,
  trace = 1e-06
}