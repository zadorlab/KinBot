    Barrier       {barrier} {reactant} {dummy}
    ! Termolecular reaction
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
    Bimolecular {product}
    ! Termolecular product
    Dummy
    End ! end Termolecular Product
