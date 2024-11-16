def start_motif(motif, natom, bond, atom, allover, eqv):
    """
    Initialize the motif search.
    If allover is >0, that atom is used as a starting point.
    """
    visit = [0] * natom
    chain = [-999] * natom
    nsteps = -1
    motifset = []
    mask = []  # to mask new patterns that allow eqv atoms in some cases
    find_motif(motif, visit, chain, nsteps, 0, -1,
               motifset, allover, natom, bond, atom, eqv, mask, True)
    return motifset


def find_motif(motif, visit, chain, nsteps, current,
               previous, motifset, allover, natom, bond, atom, eqv, mask, mask_current):
    """
    This recursive function finds a specific motif in the structure.
    FIXIT - the comments here
    nsteps: the number of steps made
    current: atom selected for testing, this is a candidate to
    be the steps'th element of motif
    previous: atom before current
    motif: a given order of given atom types connected to each other covalently
    visit: 1 if an atom was visited
    chain: atoms in order of search, when retracting, elements are overwritten.
    motifset: 2D array of the found motives
    allover < 0: all-over path-finding mode
    allover = > 0: only atom 'allover' is used for the search
    eqv: array for equvivalent atoms
    mask: mask irregulara patterns that contain equivalent atoms
    """
    if nsteps == -1:
        nsteps = 0
        previous = -1
        if allover < 0:
            for i in range(natom):
                current = i
                visit = [0] * natom
                chain = [-999] * natom
                mask_current = True
                find_motif(motif, visit, chain, nsteps, current, previous,
                           motifset, allover, natom, bond, atom, eqv, mask, mask_current)
        else:
            current = allover
            visit = [0] * natom
            find_motif(motif, visit, chain, nsteps, current, previous,
                       motifset, allover, natom, bond, atom, eqv, mask, mask_current)

    if nsteps > -1:
        if nsteps > natom:
            return 0
        # check if one of the motifs already had an equivalent atom in
        # the same position (that is not the same atom)
        eqv_list = []
        for mylist in eqv:
            if current in mylist:
                eqv_list = mylist[:]
                eqv_list.remove(current)
        for mi, m in enumerate(motifset):
            if mask[mi] == False:
                continue
            if m[nsteps] in eqv_list:
                if any([True for i in eqv_list if i in chain[:nsteps]]):
                    mask_current = False
                    break
                return 0
        if visit[current] == 1:
            return 0
        if nsteps == len(motif):
            return 0
        if current in chain[:nsteps-1]:
            return 0
        if motif[nsteps] != 'X' and atom[current] != motif[nsteps]:
            return 0
        if bond[current][previous] == 0 and previous > -1:
            return 0

        if nsteps == len(motif) - 1:
            chain[nsteps] = current
            if -999 in chain:
                motifset.append(chain[:chain.index(-999)])
            else:
                motifset.append(chain[:])
            mask.append(mask_current)
            return 0

        visit[current] = 1
        chain[nsteps] = current
        previous = current
        nsteps += 1
        for i in range(natom):
            current = i
            find_motif(motif, visit, chain, nsteps, current, previous,
                       motifset, allover, natom, bond, atom, eqv, mask, mask_current)
            visit[current] = 0

        if nsteps > 0:
            nsteps = nsteps - 1
            previous = chain[nsteps - 1]
        else:
            previous = -1
    return 0


def bondfilter(motif, bond, bondpattern, atleast=False):
    """
    For a given linear sequence of atoms it tests whether
    the bond orders match pattern bondpattern.
    E.g., bondpattern can be 1, 1, 2, 1 for a 5-long motif.
    X means any bond pattern is accepted
    atleast: the bond multiplicity is at least as in bondpattern
    E.g., the above can be 2, 1, 3, 1 as well and will be accepted
    """
    for atomi in range(len(motif)-1):
        if bondpattern[atomi] == 'X':
            continue
        if not atleast:
            if bond[motif[atomi]][motif[atomi+1]] != bondpattern[atomi]:
                return -1
        else:
            if bond[motif[atomi]][motif[atomi+1]] < bondpattern[atomi]:
                return -1
    return 0
