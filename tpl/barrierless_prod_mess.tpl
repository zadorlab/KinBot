!****************************************
    Bimolecular     {bimol_n} ! {products} 

    Fragment       {fr_n} ! {chemid fragment 1}
    RRHO
        Geometry[angstrom]            {n_atoms frag 1}
        {geom fragment 1}

        Core   RigidRotor
            SymmetryFactor            1.0
        End
    
        Frequencies[1/cm]             {n_freq1}
            {freq fragment 1}


        ElectronicLevels[1/cm]      1
            0    3

        ZeroEnergy[kcal/mol]        0.0
    End ! end fragment

    Fragment       {fr_n} ! {chemid fragment 2}
    RRHO
        Geometry[angstrom]            {n_atoms2}
        {geom fragment 2}

        Core   RigidRotor
            SymmetryFactor            1.0
        End
    
        Frequencies[1/cm]             {n_freq2}
        {freqs fragment 2}

        ElectronicLevels[1/cm]      1
            0    2

        ZeroEnergy[kcal/mol]        0.0
    End ! end fragment

        GroundEnergy[kcal/mol]      {prod_energy}

    
    End ! end bimolecular
!****************************************
