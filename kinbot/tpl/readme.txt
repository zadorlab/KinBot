mess_barrier
************

    Barrier       {rxn_name} {chemid_reac} {chemid_prod} ! {long_rxn_name}
{model}  ==> can be mess_rrho or mess_variational

mess_rrho
*********
    RRHO
        Geometry[angstrom]            {natom}
        {geom}
{core}   ==> can be mess_core_rr or mess_phasespacetheory
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

mess_phasespacetheory
**********************
    Core   PhaseSpaceTheory
        FragmentGeometry    {natom1}
        {geom1}
        FragmentGeometry    {natom2}
        {geom2}
        SymmetryFactor  {symm}
        PotentialPrefactor[au] {prefact}
        PotentialPowerExponent {exponent}
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

mess_variational
****************
    Variational
{outerts} ==> mess_2tst
{variationalmodel} ==> can be mess_rrho (repeated) w/o tunneling
{tunneling} ==> mess_tunneling

mess_2ts
********
        2TSMethod   statistical
{outerts} ==> can be mess_rrho (repeated if needed, but no tunneling) or mess_phasespacetheory

