! KinBot standalone MESS partition-function input.
TemperatureList[K]                     {temperatures}
RelativeTemperatureIncrement           0.001
AtomDistanceMin[angstrom]              0.32
Species {name}
  RRHO
    Geometry[angstrom]                 {natom}
{geom}
    Core RigidRotor
      SymmetryFactor                   {symm}
{rotconst}
    End
    Frequencies[1/cm]                  {nfreq}
{freq}
{hinderedrotor}
    ElectronicLevels[1/cm]             1
      0.0 {mult}
    ZeroEnergy[kcal/mol]               0.0
  End
End
