    Barrier       {barrier} {reactant} {dummy}
    ! Barrierless Reaction
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
 
   Bimolecular     {chemids}

{fragments}

        GroundEnergy[kcal/mol]        {ground_energy.2f}

    
    End ! end Barrierless
