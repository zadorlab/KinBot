import numpy as np
from kinbot import constants

def atom_type_table(element, nconnect, ndouble, ntriple):
    #print(f"{element}, {nconnect}, {ndouble}, {ntriple}")
    #This table does not take into account oxidation states. Neutral is assumed, unless specified in example.
    #Charges: (+)(-)

    # TODO only send lin, tri, etc., not the atom type
    atom_type = 'H'
    if element == 'H':
        if nconnect == 1:
            atom_type = 'H_lin'
    elif element == 'C':
        match nconnect:
            case 0:
                atom_type = 'C'
            case 1:
                atom_type = 'C_lin' #Ex: H - C ...
            case 2:
                match ndouble:
                    case 0:
                        atom_type = 'C_tri' #Ex: H-C:-H -> H3C.
                    case 1:
                        atom_type = 'C_lin' #Ex: O = C: -> O=C.-H
            case 3:
                match ntriple:
                    case 1:
                        atom_type = 'C_lin' #Ex: N=_C. -> N=_C-H
                    case 0:
                        match ndouble:
                            case 1:
                                atom_type = 'C_tri' #Ex: O=C.-H -> O=CH2
                            case 0:
                                atom_type = 'C_quad' #Ex: H3C. -> CH4
            case 4:  # Weird case, full valency carbon, to be double checked
                match ntriple:
                    case 1:
                        atom_type = "C_tri"
                    case 0:
                        match ndouble:
                            case 2:
                                atom_type = "C_tri"
                            case 1:
                                atom_type = "C_quad"
                            case 0:
                                atom_type = "C_tri"
    elif element == 'N':
        match nconnect:
            case 0:
                atom_type = 'N'
            case 1:
                atom_type = 'N_tri' #Ex: H - N: -> H2N.
            case 2:
                match ndouble:
                    case 0:
                        atom_type = 'N_pyr' #Ex: H2N. -> NH3
                    case 1:
                        atom_type = 'N_tri' #Ex: C=N. -> C=NH
            case 3:  # Weird case with full valency N, to be checked
                match ntriple:
                    case 1:
                        atom_type = "N_tri"
                    case 0:
                        match ndouble:
                            case 1:
                                atom_type = "N_pyr"
                            case 0:
                                atom_type = "N_quad"
            case 4:  # Weird case, more than full valency N: N+ is assumed, which is equivalent to C 
                match ntriple:
                    case 1:
                        atom_type = "N_tri"
                    case 0:
                        match ndouble:
                            case 2:
                                atom_type = "N_tri"
                            case 1:
                                atom_type = "N_quad"
                            case 0:
                                atom_type = "N_tri"
    elif element == 'O':
        match nconnect:
            case 0:
                atom_type = 'O'
            case 1:
                atom_type = 'O_tri' #Ex: H-O. -> H2O
            case 2:
                match ndouble:
                    case 0:
                        atom_type = 'O_quad' #Ex: This is weird, unless H2O(+). + R(-).
                    case 1:
                        atom_type = 'O_tri'
    elif element == 'S':
        match nconnect:
            case 0:
                atom_type = 'S'
            case 1:
                atom_type = 'S_tri' #Ex: F-S -> F2S.
            case 2:
                match ndouble:
                    case 0:
                        atom_type = 'S_pyr' #Ex: F2S -> F3S
                    case 1:
                        atom_type = 'S_lin' #Ex: S=O -> ? Not sure if the case exists.
            case 3:
                match ntriple:
                    case 1:
                        atom_type = 'S_pyr' #Ex: NS -> NSF
                    case 0:
                        match ndouble:
                            case 1:
                                atom_type = 'S_pyr' #Ex: O=S-F -> SOF2
                            case 0:
                                atom_type = 'S_bip_tri_l' #Ex: SF3 -> SF4
            case 4:
                match ntriple:
                    case 1:
                        atom_type = 'S_tri' #Not sure if this case exists. N=_S-F -> N=_SF2
                    case 0:
                        match ndouble:
                            case 2:
                                atom_type = 'S_pyr' #Ex: O=S=O -> O=SF=O
                            case 1:
                                atom_type = 'S_bip_tri_l' #Ex: O=SF2 -> O=SF3
                            case 0:
                                atom_type = 'S_bip_tri' #Ex: SF4 -> SF5
            case 5:
                match ntriple:
                    case 1:
                        match ndouble:
                            case 1:
                                atom_type = 'S_pyr'  # Ex: NSO -> NSOF
                            case 0:
                                atom_type = 'S_quad'  # Ex: NSF2 -> NSF3
                    case 0:
                        match ndouble:
                            case 2:
                                atom_type = 'S_pyr'  # Ex: SO2F -> SO2F2
                            case 1:
                                atom_type = 'S_bip_tri_l'  # Ex: SOF3 -> SOF4
                            case 0:
                                atom_type = 'S_bip_quad_t'  # Ex: SF5 (+)(-)? -> SF6
    return atom_type


def pp_length_table() -> dict[str, list[float]]:
    """Create a default list of distances to try from a given pivot point.

    Args:
        element (str): chemical symbol

    Returns:
        list[float]: list of distances for the pivot point.
    """
    pp_table: dict[str, list[float]] = {
        'H': (np.array([0.2, 0.5])*constants.BOHRtoANGSTROM).tolist(),
        'C': (np.array([0.5, 1.0])*constants.BOHRtoANGSTROM).tolist(),
        'N': (np.array([0.5, 1.0])*constants.BOHRtoANGSTROM).tolist(),
        'O': (np.array([0.5, 1.0])*constants.BOHRtoANGSTROM).tolist(),
        'S': (np.array([0.5, 1.0])*constants.BOHRtoANGSTROM).tolist(),
        'X': (np.array([0.5, 1.0])*constants.BOHRtoANGSTROM).tolist()}

    return pp_table
