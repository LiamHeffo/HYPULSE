-- Driver CEA Gas model made automatically by PITOT3 on the 15/12/2023 01:48:40

model = "CEAGas"

CEAGas = {
  mixtureName = "driver_gas",
  speciesList = {'O2', 'Ar', 'H2'},
  reactants = {O2 = 18.3, H2 = 36.7, Ar = 45},
  inputUnits = 'moles',
  withIons = false,
  trace = 1e-06
}