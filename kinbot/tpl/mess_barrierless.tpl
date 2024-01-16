  Bimolecular  {chemids}
{fragments}
    GroundEnergy[kcal/mol]        {ground_energy}
    End ! Barrierless

  Barrier       {barrier} {reactant} {prod} 
    RRHO
      Stoichiometry {stoich}
        Core PhaseSpaceTheory
{frag1}
{frag2}
          SymmetryFactor         {pstsymm}
          PotentialPrefactor[au] 0.6
          PotentialPowerExponent 6. 
        End ! Core
      Frequencies[1/cm]             {nfreq}
{freq}
{hinderedrotor}
        ElectronicLevels[1/cm]      1
            0.    {mult} 
        ZeroEnergy[kcal/mol]      {ground_energy} 
    End ! RRHO

!****************************************

