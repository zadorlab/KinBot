mess_barrier
************

    Barrier       {rxn_name} {chemid_reac} {chemid_prod} ! {long_rxn_name}
{model}  ==> mess_rrho

mess_rrho
*********
    RRHO
        Geometry[angstrom]            {natom}
        {geom}
{core}   ==> mess_core_rr
        Frequencies[1/cm]             {nfreq}
            {freq}
{rotor}   ==> mess_hinderedrotor, which is now just the 1D model, later can be multirotor as well
{tunneling} ==> mess_tunneling, which is just the Eckart model, but should be fine
        ElectronicLevels[1/cm]      {nelec}
            {charge}    {mult}
        ZeroEnergy[kcal/mol]        {zeroenergy}
    End ! RRHO


mess_core_rr
************
        Core   RigidRotor
            SymmetryFactor            {symm}
        End ! Core

mess_hindered_rotor
*******************
        Rotor     Hindered
        ThermalPowerMax         50.
            Group                     {group}
            Axis                      {axis}
            Symmetry                  {rotorsymm}
            Potential[kcal/mol]       {nrotorpot}
                {rotorpot}
        End ! Rotor

mess_tunneling
**************
        Tunneling   Eckart
                CutoffEnergy[kcal/mol]    {cutoff}
                ImaginaryFrequency[1/cm]  {imfreq}
                WellDepth[kcal/mol]       {welldepth1}
                WellDepth[kcal/mol]       {welldepth2}
        End ! Tunneling

