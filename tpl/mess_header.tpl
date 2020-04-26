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
      Epsilons[1/cm]                    {e_coll} {e_well.3f}
      Sigmas[angstrom]                  {s_coll} {s_well}
      Masses[amu]                       {m_coll} {m_well}
    End
