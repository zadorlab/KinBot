import os
import subprocess
import time
import xml.etree.cElementTree as ET
import xml.dom.minidom as minidom

from kinbot import kb_path
from kinbot import constants


class MESMER:
    """
    Class that reads and write all MESMER files and runs MESMER
    """
    def __init__(self, par, species):
        self.par = par
        self.species = species

    def write_input(self):
        """
        write the input for all the wells, bimolecular products and barriers
        both in a separate file, as well as in one large ME file
        """
        # filter ts's with the same reactants and products:
        ts_unique = {}  # key: ts name, value: [prod_name, energy]
        for index, reaction in enumerate(self.species.reac_obj):
            if self.species.reac_ts_done[index] == -1:
                prod_name = '_'.join([str(pi.chemid) for pi in reaction.products])
                energy = reaction.ts.energy
                new = 1
                remove = []
                for ts in ts_unique:
                    if ts_unique[ts][0] == prod_name:
                        # check for the barrier with the lowest energy
                        if ts_unique[ts][1] > energy:
                            # remove the current barrier
                            remove.append(ts)
                        else:
                            new = 0
                for ts in remove:
                    ts_unique.pop(ts, None)
                if new:
                    ts_unique[reaction.instance_name] = [prod_name, energy]

        root = ET.Element('me:mesmer',
                          {'xmlns': 'http://www.xml-cml.org/schema',
                           'xmlns:me': 'http://www.chem.leeds.ac.uk/mesmer',
                           'xmlns:xsi': 'http://www.w3.org/2001/XMLSchema-instance'})

        ET.SubElement(root, 'me:title').text = str(self.species.chemid)
        mollist = ET.SubElement(root, 'moleculeList')
        reaclist = ET.SubElement(root, 'reactionList')

        # write the mess input for the different blocks
        wells = [self.species.chemid]

        # write the initial reactant
        well0_energy = self.species.energy - self.species.zpe
        self.write_well(self.species, mollist, well0_energy)

        for index, reaction in enumerate(self.species.reac_obj):
            if reaction.instance_name in ts_unique:
                # get the energy of the products compared to the well0
                prod_energy = 0.0 - well0_energy
                # get the produt with the most atoms and assign the energy to that species
                # the other get an energy of 0
                max_natom = 0
                max_natom_chemid = 0
                """
                for st_pt_opt in reaction.prod_opt:
                    prod_energy += st_pt_opt.species.energy
                    prod_energy += st_pt_opt.species.zpe
                    if st_pt_opt.species.natom > max_natom:
                        max_natom = st_pt_opt.species.natom
                        max_natom_chemid = st_pt_opt.species.chemid
                """
                for st_pt_opt in reaction.prod_opt:
                    st_pt = st_pt_opt.species
                    prod_energy = st_pt_opt.species.energy
                    prod_energy += st_pt_opt.species.zpe
                    if st_pt.chemid not in wells:
                        energy = 0.
                        if st_pt.chemid == max_natom_chemid:
                            energy = prod_energy
                        self.write_well(st_pt, mollist, prod_energy)
                        wells.append(st_pt.chemid)
                self.write_barrier(reaction, mollist, reaclist)

        # add a bath gas molecule
        molecule = ET.SubElement(mollist, 'molecule', {'id': 'He'})
        propertylist = ET.SubElement(molecule, 'propertyList')
        # add epsilon
        epsilon = ET.SubElement(propertylist, 'property', {'dictRef': 'me:epsilon'})
        ET.SubElement(epsilon, 'scalar').text = '48.0'
        # add sigma
        sigma = ET.SubElement(propertylist, 'property', {'dictRef': 'me:sigma'})
        ET.SubElement(sigma, 'scalar').text = '3.90'
        # add the molecular weight
        mw = ET.SubElement(propertylist, 'property', {'dictRef': 'me:mw'})
        ET.SubElement(mw, 'scalar').text = '28.0'

        # add the conditions
        conditions = ET.SubElement(root, 'me:conditions')
        ET.SubElement(conditions, 'me:bathGas').text = 'He'
        pts = ET.SubElement(conditions, 'me:PTs')
        for p in self.par['PressureList']:
            for t in self.par['TemperatureList']:
                ET.SubElement(pts,
                              'me:PTpair',
                              {'units': 'Torr',
                               'p': str(p),
                               't': str(t),
                               'precision': 'd'})

        # add the model parameters
        modelparams = ET.SubElement(root, 'me:modelParameters')
        ET.SubElement(modelparams, 'me:grainSize', {'units': 'cm-1'}).text = '50'
        ET.SubElement(modelparams, 'me:energyAboveTheTopHill').text = '30.0'

        # add the control parameters
        control = ET.SubElement(root, 'me:control')
        ET.SubElement(control, 'me:printSpeciesProfile')
        ET.SubElement(control, 'me:testRateConstants')
        ET.SubElement(control, 'me:eigenvalues').text = '3'

        st = ET.tostring(root, 'utf-8')
        st = minidom.parseString(st)
        with open('me/mesmer.xml', 'w') as fout:
            fout.write(st.toprettyxml(indent=' '))

        return 0

    def write_well(self, species, mollist, energy):
        """
        Create the block for MESS for a well.
        """
        molecule = ET.SubElement(mollist, 'molecule', {'id': str(species.chemid)})
        atomarray = ET.SubElement(molecule, 'atomArray')
        for i, at in enumerate(species.atom):
            args = {'id': 'a{}'.format(i+1)}
            args['elementType'] = at
            args['x3'] = '{:.8f}'.format(species.geom[i][0])
            args['y3'] = '{:.8f}'.format(species.geom[i][1])
            args['z3'] = '{:.8f}'.format(species.geom[i][2])
            at = ET.SubElement(atomarray, 'atom', args)

        bondarray = ET.SubElement(molecule, 'bondArray')
        bond_id = 1
        bond_ref = {}
        for i in range(species.natom-1):
            for j in range(i+1, species.natom):
                if species.bond[i][j] > 0:
                    bond_name = 'b{}'.format(bond_id)
                    bond_ref['{}_{}'.format(i, j)] = bond_name
                    args = {'id': 'b{}'.format(bond_id)}
                    args['atomRefs2'] = "a{} a{}".format(i+1, j+1)
                    args['order'] = "{}".format(species.bond[i][j])
                    b = ET.SubElement(bondarray, 'bond', args)
                    bond_id += 1
        propertylist = ET.SubElement(molecule, 'propertyList')

        # add the zpe
        zpe = ET.SubElement(propertylist, 'property', {'dictRef': 'me:ZPE'})
        ET.SubElement(zpe, 'scalar', {'units': 'Hartree'}).text = str(energy)

        # add the multiplicity
        mult = ET.SubElement(propertylist, 'property', {'dictRef': 'me:spinMultiplicity'})
        ET.SubElement(mult, 'scalar').text = str(species.mult)

        # add the external symmetry number
        sigma = ET.SubElement(propertylist, 'property', {'dictRef': 'me:symmetryNumber'})
        ET.SubElement(sigma, 'scalar').text = str(species.sigma_ext / species.nopt)
        
        # add the vibrational frequencies
        if len(species.reduced_freqs) > 0:
            vibs = ET.SubElement(propertylist, 'property', {'dictRef': 'me:vibFreqs'})
            ET.SubElement(vibs, 'array', {'units': 'cm-1'}).text = ' '.join([str(fi) for fi in species.reduced_freqs])

            # add the rotor potentials
            if self.par['rotor_scan'] and len(species.dihed) > 0:
                qmrotors = ET.SubElement(molecule, 'me:DOSCMethod', {'xsi:type': 'me:QMRotors'})
                for i, rot in enumerate(species.dihed):
                    hinderedrotor = ET.SubElement(molecule, 'me:ExtraDOSCMethod', {'xsi:type': 'me:HinderedRotorQM1D'})
                    ET.SubElement(hinderedrotor, 'me:bondRef').text = bond_ref['_'.join(sorted([str(rot[1]), str(rot[2])]))]
                    pot = ET.SubElement(hinderedrotor,
                                        'me:HinderedRotorPotential',
                                        {'format': 'numerical',
                                         'units': 'kcal/mol',
                                         'expansionSize': '7',
                                         'useSineTerms': 'yes'})
                    ens = species.hir.hir_energies[i]
                    rotorpot = [(ei - ens[0])*constants.AUtoKCAL for ei in ens]
                    for j, ri in enumerate(rotorpot):
                        ET.SubElement(pot,
                                      'me:PotentialPoint',
                                      {'angle': str(360/species.hir.nrotation*j),
                                       'potential': str(ri)})
                    ET.SubElement(pot,
                                  'me:PotentialPoint',
                                  {'angle': str(360.0),
                                   'potential': str(rotorpot[0])})
                    ET.SubElement(hinderedrotor, 'me:periodicity').text = str(species.sigma_int[rot[1]][rot[2]])
            else:
                qmrotors = ET.SubElement(molecule, 'me:DOSCMethod', {'name': 'QMRotors'})

        # add the energy transfer model
        etransfer = ET.SubElement(molecule, 'me:energyTransferModel', {'xsi:type': 'me:ExponentialDown'})
        ET.SubElement(etransfer, 'me:deltaEDown', {'units': 'cm-1'}).text = '-1'

    def write_barrier(self, reaction, mollist, reaclist):
        """
        Create the block for a MESS barrier.
        """
        rxn = ET.SubElement(reaclist, 'reaction', {'id': reaction.instance_name})

        # add the reactant refs
        reactant = ET.SubElement(rxn, 'reactant')
        ET.SubElement(reactant, 'molecule', {'ref': str(self.species.chemid), 'role': 'modelled'})

        # add the product refs
        if len(reaction.products) == 1:
            product = ET.SubElement(rxn, 'product')
            ET.SubElement(product, 'molecule', {'ref': str(reaction.products[0].chemid), 'role': 'modelled'})
        else:
            for st_pt in reaction.products:
                product = ET.SubElement(rxn, 'product')
                ET.SubElement(product, 'molecule', {'ref': str(st_pt.chemid), 'role': 'sink'})

        # add the transition state ref
        ts = ET.SubElement(rxn, 'me:transitionState')
        ET.SubElement(ts, 'molecule', {'ref': reaction.instance_name, 'role': 'transitionState'})

        # add the MCRCMethod
        ET.SubElement(rxn, 'me:MCRCMethod', {'name': 'RRKM'})

        # add the tunneling
        barriers = [
            ((reaction.ts.energy + reaction.ts.zpe) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL,
            ((reaction.ts.energy + reaction.ts.zpe) - sum([(opt.species.energy + opt.species.zpe) for opt in reaction.prod_opt])) * constants.AUtoKCAL,
        ]
        if all([bi > 0 for bi in barriers]):
            ET.SubElement(rxn, 'me:tunneling', {'name': 'Eckart'})

        # add the transition state to the molecule list
        molecule = ET.SubElement(mollist, 'molecule', {'id': reaction.instance_name})
        atomarray = ET.SubElement(molecule, 'atomArray')
        for i, at in enumerate(reaction.ts.atom):
            args = {'id': 'a{}'.format(i+1)}
            args['elementType'] = at
            args['x3'] = '{:.8f}'.format(reaction.ts.geom[i][0])
            args['y3'] = '{:.8f}'.format(reaction.ts.geom[i][1])
            args['z3'] = '{:.8f}'.format(reaction.ts.geom[i][2])
            at = ET.SubElement(atomarray, 'atom', args)

        bondarray = ET.SubElement(molecule, 'bondArray')
        bond_id = 1
        bond_ref = {}
        for i in range(reaction.ts.natom-1):
            for j in range(i+1, reaction.ts.natom):
                if reaction.ts.bond[i][j] > 0:
                    bond_name = 'b{}'.format(bond_id)
                    bond_ref['{}_{}'.format(i, j)] = bond_name
                    args = {'id': 'b{}'.format(bond_id)}
                    args['atomRefs2'] = "a{} a{}".format(i+1, j+1)
                    args['order'] = "{}".format(reaction.ts.bond[i][j])
                    b = ET.SubElement(bondarray, 'bond', args)
                    bond_id += 1
        propertylist = ET.SubElement(molecule, 'propertyList')

        # add the zpe
        if self.par['pes']:
            energy = '{zeroenergy}'
        else:
            energy = reaction.ts.energy + reaction.ts.zpe
        zpe = ET.SubElement(propertylist, 'property', {'dictRef': 'me:ZPE'})
        ET.SubElement(zpe, 'scalar', {'units': 'Hartree'}).text = str(energy)

        # add the multiplicity
        mult = ET.SubElement(propertylist, 'property', {'dictRef': 'me:spinMultiplicity'})
        ET.SubElement(mult, 'scalar').text = str(reaction.ts.mult)

        # add the external symmetry number
        sigma = ET.SubElement(propertylist, 'property', {'dictRef': 'me:symmetryNumber'})
        ET.SubElement(sigma, 'scalar').text = str(reaction.ts.sigma_ext / reaction.ts.nopt)

        # add the vibrational frequencies
        vibs = ET.SubElement(propertylist, 'property', {'dictRef': 'me:vibFreqs'})
        ET.SubElement(vibs, 'array', {'units': 'cm-1'}).text = ' '.join([str(fi) for fi in reaction.ts.reduced_freqs[1:]])

        # add the imaginary frequencies
        imfreq = ET.SubElement(propertylist, 'property', {'dictRef': 'me:imFreqs'})
        ET.SubElement(imfreq, 'scalar', {'units': 'cm-1'}).text = str(-reaction.ts.reduced_freqs[0])

        # add the rotor potentials
        if self.par['rotor_scan'] and len(reaction.ts.dihed) > 0:
            qmrotors = ET.SubElement(molecule, 'me:DOSCMethod', {'xsi:type': 'me:QMRotors'})
            for i, rot in enumerate(reaction.ts.dihed):
                hinderedrotor = ET.SubElement(molecule, 'me:ExtraDOSCMethod', {'xsi:type': 'me:HinderedRotorQM1D'})
                ET.SubElement(hinderedrotor, 'me:bondRef').text = bond_ref['_'.join(sorted([str(rot[1]), str(rot[2])]))]
                pot = ET.SubElement(hinderedrotor, 'me:HinderedRotorPotential', {'format': 'numerical', 'units': 'kcal/mol', 'expansionSize': '7', 'useSineTerms': 'yes'})
                ens = reaction.ts.hir.hir_energies[i]
                rotorpot = [(ei - ens[0])*constants.AUtoKCAL for ei in ens]
                for j, ri in enumerate(rotorpot):
                    ET.SubElement(pot, 'me:PotentialPoint', {'angle': str(360/reaction.ts.hir.nrotation*j), 'potential': str(ri)})
                ET.SubElement(pot, 'me:PotentialPoint', {'angle': str(360.0), 'potential': str(rotorpot[0])})
                ET.SubElement(hinderedrotor, 'me:periodicity').text = str(reaction.ts.sigma_int[rot[1]][rot[2]])
        else:
            qmrotors = ET.SubElement(molecule, 'me:DOSCMethod', {'name': 'QMRotors'})

        # add the energy transfer model
        etransfer = ET.SubElement(molecule, 'me:energyTransferModel', {'xsi:type': 'me:ExponentialDown'})
        ET.SubElement(etransfer, 'me:deltaEDown', {'units': 'cm-1'}).text = '-1'

    def run(self):
        """
        write a pbs/slurm file for the me/all.xml input file
        submit the pbs/slurm file to the queue
        wait for the mess run to finish
        """

        # open the the header and the specific templates
        if self.par['queue_template'] == '':
            q_file = f'{kb_path}/tpl/{self.par["queuing"]}.tpl'
        else:
            q_file = self.par['queue_template']
        with open(q_file) as f:
            tpl_head = f.read()
        q_file = f'{kb_path}/tpl/{self.par["queuing"]}_mesmer.tpl'
        with open(q_file) as f:
            tpl = f.read()
        submitscript = 'run_mesmer' + constants.qext[self.par['queuing']] 
        with open(submitscript, 'a') as qu:
            if self.par['queue_template'] == '':
                if self.par['queuing'] == 'pbs':
                    qu.write((tpl_head + tpl).format(name='mesmer', ppn=self.par['ppn'], queue_name=self.par['queue_name'], dir='me'))
                elif self.par['queuing'] == 'slurm':
                    qu.write((tpl_head + tpl).format(name='mesmer', ppn=self.par['ppn'], queue_name=self.par['queue_name'], dir='me'), slurm_feature='')
            else:
                qu.write(tpl_head)
                qu.write(tpl)

        command = [constants.qsubmit[self.par['queuing']], submitscript]
        process = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        out = out.decode()
        if self.par['queuing'] == 'pbs':
            pid = out.split('\n')[0].split('.')[0]
        elif self.par['queuing'] == 'slurm':
            pid = out.split('\n')[0].split()[-1]

        while 1:  
            devnull = open(os.devnull, 'w')
            if self.par['queuing'] == 'pbs':
                command = 'qstat -f | grep ' + '"Job Id: ' + pid + '"' + ' > /dev/null'
            elif self.par['queuing'] == 'slurm':
                command = 'scontrol show job ' + pid + ' | grep "JobId=' + pid + '"' + ' > /dev/null'
            if int(subprocess.call(command, shell=True, stdout=devnull, stderr=devnull)) == 0:
                time.sleep(1)
            else:
                break
        return 0



