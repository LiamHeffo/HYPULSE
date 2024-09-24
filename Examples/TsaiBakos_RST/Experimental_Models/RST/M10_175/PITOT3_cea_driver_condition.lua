-- Driver CEA Gas model made automatically by PITOT3 on the 15/12/2023 01:51:27

model = "CEAGas"

CEAGas = {
  mixtureName = "driver_gas",
  speciesList = {'O2', 'Ar', 'H2'},
  reactants = {O2 = 23.3, H2 = 46.7, Ar = 30},
  inputUnits = 'moles',
  withIons = false,
  trace = 1e-06
}