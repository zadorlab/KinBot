from fragments.py import Fragment

class VRC_TST_surfaces()
'''
Class that generates a list of pivot points and a list of pivot points
and their associated distances depending on the intermolecular distance.
'''

    def __init__(self, distances=None, surfaces=None, nfrag=None)
        self.distances = []
        self.surfaces = []
        self.nfrag = len(fragments)

    def create_surfaces(self, par, fragments)
        self.distances.append(i for i in arange(par['barrierless_saddle_start'],
                                                par['barrierless_saddle_stop'],
                                                par['barrierless_saddle_step']))
        for this_frag in fragments:
            #TODO:function of the fragment class to be written
            #Must return the number of the atom in the fragment correponding to the number in the reactant.
            this_frag.set_labels(frag_number, self.nfrag)

        #Recover the index of the reactive atoms from the reaction name.
        summary = open(this_frag.parent_chemid + '/summary_' + this_frag.parent_chemid + '.out', 'r').readlines()
        #Check if the line corresponds to the reaction before taking the indices
        for line in summary:
            if line.startswith('SUCCESS') and 'hom_sci' in line:
                pieces = line.split()
                current_parent = pieces[2].split('_')[0]

                corresponds = 0
                #If this reaction has a number of product different than the number of fragments,
                #then skip this iteration.
                if len(pieces[3:]) == len(fragments):
                    pass
                else
                    continue

                for current_product, this_frag in zip(pieces[3:], fragments):
                    if current_parent == this_frag.parent_chemid and current_product == this_frag.chemid:
                        corresponds = 1
                    else
                        corresponds = 0
                        break
                #When the reaction is found, create the list of reactive atoms        
                if corresponds:
                    reactive_atoms = []
                    for i in pieces[2].split('_')[3:]:
                        reactive_atoms.append(int(i))
                    break

        #for all distances
        for dist in self.distances:
            #Create a coordinate system with all the fragments
            #TODO: write the global_system class
            assembled = global_system.create_sys(fragments)
            #Set the position of the pivot points in the fragment object
            #TODO: write the set_pivot_point() method in the class global_system
            #It must call the creation of a pivot point attribute of the fragment class
            assembled.set_pivot_points(reactive_atoms)
            #Create the "surface" attribute (dictionary) as expected by rotd_py
            #TODO Must return the coordinates of the pivot points from the fragment object,
            #and a pivot point distances matrix computed in the global_system class.
            assembled.set_surface()
            
            self.surfaces.append(assembled.get_surface())

        return self.surfaces
