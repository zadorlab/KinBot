TemperatureList[K]                      {TemperatureList}
PressureList[torr]                      {PressureList}
EnergyStepOverTemperature               {EnergyStepOverTemperature}
ExcessEnergyOverTemperature             {ExcessEnergyOverTemperature}
ModelEnergyLimit[kcal/mol]              {ModelEnergyLimit}
CalculationMethod                       {CalculationMethod}
ChemicalEigenvalueMax                   {ChemicalEigenvalueMax}
Reactant                                {Reactant}
MicroRateOutput micro.out
MicroEnerMax[kcal/mol] 100.
MicroEnerMin[kcal/mol] 0.
MicroEnerStep[kcal/mol] 1.
Model
  EnergyRelaxation
    Exponential
      Factor[1/cm]                      {EnergyRelaxationFactor}
      Power                             {EnergyRelaxationPower}
      ExponentCutoff                    {EnergyRelaxationExponentCutoff}
    End
  CollisionFrequency
    LennardJones
      Epsilons[1/cm]                    {e_coll} {e_well}
      Sigmas[angstrom]                  {s_coll} {s_well}
      Masses[amu]                       {m_coll} {m_well}
    End
