TemperatureList[K]                      {TemperatureList}
PressureList[torr]                      {PressureList}
EnergyStepOverTemperature               {EnergyStepOverTemperature}
ExcessEnergyOverTemperature             {ExcessEnergyOverTemperature}
ModelEnergyLimit[kcal/mol]              {ModelEnergyLimit}
CalculationMethod                       {CalculationMethod}
ChemicalEigenvalueMax                   {ChemicalEigenvalueMax}
Reactant                                {Reactant}
Model
  EnergyRelaxation
    Exponential
      Factor[1/cm]                      {EnergyRelaxationFactor}
      Power                             {EnergyRelaxationPower}
      ExponentCutoff                    {EnergyRelaxationExponentCutoff}
    End
  CollisionFrequency
    LennardJones
      Epsilons[1/cm]                    {Epsilons}
      Sigmas[angstrom]                  {Sigmas}
      Masses[amu]                       {Masses}
    End
