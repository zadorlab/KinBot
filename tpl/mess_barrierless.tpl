    Barrier       {barrier} {reactant} {prod} {dummy}
    RRHO
        Geometry[angstrom]            2
	    O 1.0 0.0 0.0
            O -1.0 0.0 0.0

        Core   RigidRotor
            SymmetryFactor            1
        End
    
        Frequencies[1/cm]             1
            100.00
        
        ElectronicLevels[1/cm]      1
            0    1

        ZeroEnergy[kcal/mol]        200.
    End ! end barrier

!****************************************
   Bimolecular  {prod} {dummy}     

{fragments}

        GroundEnergy[kcal/mol]        {ground_energy}

    
    End ! end Barrierless
