
def atom_type_table(element, nconnect, ndouble, ntriple):
    #This table does not take into account oxidation states. Neutral is assumed, unless specified in example.
    #Charges: (+)(-)
    if element == 'H':
        atom_type = 'H'
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
    elif element == 'O':
        match nconnect:
            case 0:
                atom_type = 'O'
            case 1:
                atom_type = 'O_tri' #Ex: H-O. -> H2O
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
                                atom_type = 'S_pyr' #Ex: NSO -> NSOF
                            case 0:
                                atom_type = 'S_quad' #Ex: NSF2 -> NSF3
                    case 0:
                        match ndouble:
                            case 2:
                                atom_type = 'S_pyr' #Ex: SO2F -> SO2F2
                            case 1:
                                atom_type = 'S_bip_tri_l' #Ex: SOF3 -> SOF4
                            case 0:
                                atom_type = 'S_bip_quad_t' #Ex: SF5 (+)(-)? -> SF6
    else :
        atom_type = 'dummy'
    return atom_type

def pp_lenght_table(element):
    match element:
        case 'H':
            return 0.5
        case 'C':
            return 0.5
        case 'N':
            return 0.5
        case 'O':
            return 0.5
        