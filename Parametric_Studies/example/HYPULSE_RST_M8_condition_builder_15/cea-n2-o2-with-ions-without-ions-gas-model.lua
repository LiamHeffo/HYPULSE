-- modified gas model made by PITOT3 on the 24/06/2024 17:00:00 without ions for low temperature operation

model = "CEAGas"

CEAGas = {
  mixtureName = "n2-o2-with-ions-without-ions",
  speciesList = {'N2', 'O2', 'NO', 'N2O', 'NO2', 'N', 'O'},
  reactants = {N2 = 0.79, O2 = 0.21},
  inputUnits = 'moles',
  withIons = false,
  trace = 1e-06
}