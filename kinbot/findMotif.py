###################################################
##                                               ##
## This file is part of the KinBot code v2.0     ##
##                                               ##
## The contents are covered by the terms of the  ##
## BSD 3-clause license included in the LICENSE  ##
## file, found at the root.                      ##
##                                               ##
## Copyright 2018 National Technology &          ##
## Engineering Solutions of Sandia, LLC (NTESS). ##
## Under the terms of Contract DE-NA0003525 with ##
## NTESS, the U.S. Government retains certain    ##
## rights to this software.                      ##
##                                               ##
## Authors:                                      ##
##   Judit Zador                                 ##
##   Ruben Van de Vijver                         ##
##                                               ##
###################################################
    def find_motif(steps, x, y, motif, bond, visit, path, finpath, nmotif, zz):
        """  
        This recursive function finds a specific motif in the structure. 
        A motif is defined as a given order of given atom types connected to each other covalently.
        steps is the number of steps we made in the current recursion
    
        x is the current atom selected for testing, this is a candidate to be the steps'th element of motif
        y is the atom at the previous step
        motif contains the atoms in the motif. 
        bond is the bond matrix
        visit holds 1 if an atom was visited, otherwise the element is 0.
        path holds the current chain of atoms (as far as we got). When retracting, elements are overwritten.
        finpath is the 2D array of the found pathways (to save and pass back to the main code) 
        nmotif is the number of motifs found
        zz < 0: it is a regular, all-over path-finding mode
        zz = > 0: it is a mode, when only one starting point is used for the search (a given atom). Its number is zz.
        """
    
    
        if steps > -1: # here were are already in the recursion part
        # conditions to step out from the loop in order
            if steps > natom: return 0 # we run out of atom (can happen when the motif is longer than the longest path in the structure) 
            if visit[x] == 1: return 0 # we've been here already. 
            if motif[steps] != 'X' and element[x] !=  motif[steps]: return 0 # atom[x] is not a wild card and is the wrong atom type
            if bond[x][y] == 0 and y > -1: return 0 # x is not connected to y, the previous atom in the search (and there was a previous) 
    
            # if none of the above conditions are true, it means that we found a good atom x, which needs to be marked 
            if motif[steps+1] == '0': # we reached the end of the motif, time to retract
                path[steps] = x
                for i in range(steps+1):
                    finpath[nmotif][i] = path[i]
                nmotif = nmotif+1
                return 0
    
            visit[x] = 1 # we mark the current atom, x, as visited in this particular path
            path[steps] = x    
            y = x # the current atom is the "previous" in the next loop
            steps = steps+1;
            for i in range(natom):
                x = i;
                find_motif ( steps, x, y, motif, bond, visit, path, finpath, nmotif, zz )
    
            if steps > 0: steps = steps-1
            if steps > 0: y = path[steps-1] # one before the actual
            else: y = -1
            visit[x] = 0
    
    
        # initialization of the recursive loop 
        if steps == -1:
            # initialize arrays
            motif = ["0" for x in range(natom)]  
            motif[0:len(mot)] = mot
            path = ["0" for x in range(natom)]
            
            steps = 0
            y = -1 # "previous atom" 
            if zz < 0:
                for i in range(natom): # this steps through each atom
                    x = i # i could be the counter as well, but it is easier to follow the logic
                    visit =  np.zeros((natom,), dtype=np.int) # initialize the visit array at the head of the graph
    
                    find_motif (steps, x, y, motif, bond, visit, path, finpath, nmotif, zz) # here we descend into the next level
    
            else:
                visit =  np.zeros((natom,), dtype=np.int) # initialize the visit array at the head of the graph
    
                find_motif (steps, zz, y, motif, bond, visit, path, finpath, nmotif, zz) # here we descend into the next level
    
        return 0
    
    
    
    
