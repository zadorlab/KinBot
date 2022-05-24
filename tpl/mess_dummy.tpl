  Barrier       {barrier} {reactant} {dummy}
    RRHO
      Geometry[angstrom]            2
      Core   RigidRotor
        SymmetryFactor            1
      End
      Frequencies[1/cm]             1
            100.00
        ElectronicLevels[1/cm]      1
            0    1
        ZeroEnergy[kcal/mol]        200.
    End ! RRHO

!****************************************
    Bimolecular {dummy}
      Dummy
      End ! Bimolecular
