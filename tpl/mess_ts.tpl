    Barrier       {rxn_name} {chemid_reac} {chemid_prod} ! {long_rxn_name}
    RRHO
        Geometry[angstrom]            {natom}
	    {geom}

        Core   RigidRotor
            SymmetryFactor            {symm}
        End
    
        Frequencies[1/cm]             {nfreq}
            {freq}
        
{hinderedrotor}

{tunneling}
        
        ElectronicLevels[1/cm]      {nelec}
            {charge}    {mult}

        ZeroEnergy[kcal/mol]        {zeroenergy.2f}
    End ! end barrier
