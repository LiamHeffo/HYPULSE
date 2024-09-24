-- Driver CEA Gas model made automatically by PITOT3 on the 20/01/2024 15:53:28

model = "CEAGas"

CEAGas = {
  mixtureName = "driver_gas",
  speciesList = {'O2', 'Ar', 'H2'},
  reactants = {O2 = 13.3, H2 = 26.7, Ar = 60},
  inputUnits = 'moles',
  withIons = false,
  trace = 1e-06
}