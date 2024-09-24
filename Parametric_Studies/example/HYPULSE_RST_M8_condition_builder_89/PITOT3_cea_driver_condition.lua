-- Driver CEA Gas model made automatically by PITOT3 on the 24/06/2024 17:12:20

model = "CEAGas"

CEAGas = {
  mixtureName = "driver_gas",
  speciesList = {'O2', 'Ar', 'H2'},
  reactants = {O2 = 13.3, H2 = 26.7, Ar = 60},
  inputUnits = 'moles',
  withIons = false,
  trace = 1e-06
}