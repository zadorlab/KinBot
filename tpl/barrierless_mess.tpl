!****************************************
Barrier         RO2_TS          w_1 {bless_prod} ! {bless_rxn}
  RRHO
    Stoichiometry               {stoichiometry}
    Core Rotd
      File                      {mc_flux.out}
      SymmetryFactor  8
      End
    Frequencies[1/cm]  {n_freq}
        {freqs fragment 1}
        {freqs fragment 2}

    ! no tunneling bc no barrier

    ZeroEnergy[kcal/mol]                {prod_energy}
    ! ElectronicEnergy[kcal/mol]          0
    ElectronicLevels[1/cm]              1
          0     2
End ! end barrierless reaction
!****************************************
