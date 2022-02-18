import os
import numpy as np
import copy
import time
import pkg_resources
from kinbot import modify_geom
from kinbot import geometry
from reactions.reac_abstraction import abstraction_align


def carry_out_reaction(rxn, step, command, bimol=0):
    """
    Verify what has been done and what needs to be done
    skip: boolean which tells to skip the first 12 steps in case of an instance shorter than 4
    scan: boolean which tells if this is part of an energy scan along a bond length coordinate
    """
    if step > 0:
        status = rxn.qc.check_qc(rxn.instance_name)
        if status != 'normal' and status != 'error': return step
  
    kwargs = rxn.qc.get_qc_arguments(rxn.instance_name, rxn.species.mult, rxn.species.charge, ts=1,
                                     step=step, max_step=rxn.max_step, scan=rxn.scan)
    if step == 0:
        if rxn.qc.is_in_database(rxn.instance_name):
            if rxn.qc.check_qc(rxn.instance_name) == 'normal': 
                err, freq = rxn.qc.get_qc_freq(rxn.instance_name, rxn.species.natom)
                if err == 0 and len(freq) > 0.:
                    err, geom = rxn.qc.get_qc_geom(rxn.instance_name, rxn.species.natom)
                    step = rxn.max_step + 1
                    return step
        if rxn.skip and len(rxn.instance) < 4:
            step = 12
        geom = rxn.species.geom
        if bimol:
            if rxn.family_name == 'abstraction':
                geom = abstraction_align(rxn.species.geom, rxn.instance, rxn.species.fragA.natom)
#                # align O----H, O is at origin, H is on +Z axis
#                g0 = copy.deepcopy(geometry.translate_and_rotate(rxn.species.geom, rxn.instance[0], rxn.instance[1]))
#                H_pos = g0[rxn.instance[1]]  # position of the H atom
#                # align H-C, H is at origin, C is on +Z axis at an angle
#                g1 = copy.deepcopy(geometry.translate_and_rotate(rxn.species.geom, rxn.instance[1], rxn.instance[2], setangle=20.))
#                g1 += H_pos  # shift the whole thing
#                geom = np.concatenate((g0[:rxn.species.fragA.natom], g1[rxn.species.fragA.natom:]))  # does not work for reverse
#                for bondmx in rxn.species.bonds:  # add missing bond to detect the need for dummy
#                    bondmx[rxn.instance[0]][rxn.instance[1]] = 1
#                    bondmx[rxn.instance[1]][rxn.instance[0]] = 1

    elif step == rxn.max_step and rxn.scan:
        err, geom = rxn.qc.get_qc_geom(rxn.instance_name, rxn.species.natom, allow_error=1, previous=1)
    else:
        err, geom = rxn.qc.get_qc_geom(rxn.instance_name, rxn.species.natom, allow_error=1)

    step, fix, change, release = rxn.get_constraints(step, geom)

    if step > rxn.max_step:
        return step
    
    #apply the geometry changes here and fix the coordinates that changed
    change_starting_zero = []
    for c in change:
        c_new = [ci - 1 for ci in c[:-1]]
        c_new.append(c[-1])
        change_starting_zero.append(c_new)
    if len(change_starting_zero) > 0:
        success, geom = modify_geom.modify_coordinates(rxn.species, rxn.instance_name, geom, change_starting_zero, rxn.species.bond)
        for c in change:
            fix.append(c[:-1])
        change = []

    #atom, geom, dummy = rxn.qc.add_dummy(rxn.species.atom, geom, rxn.species.bond)

    kwargs['addsec'] = ''
    for fixi in fix:
        kwargs['addsec'] += f"{' '.join(str(f) for f in fixi)} F\n"
    for chi in change:
        kwargs['addsec'] += f"{' '.join(str(ch) for ch in changei)} F\n"
    for reli in release:
        kwargs['addsec'] += f"{' '.join(str(rel) for rel in reli)} A\n"

    if step < rxn.max_step:
        template_file = pkg_resources.resource_filename('tpl', 'ase_{qc}_ts_search.tpl.py'.format(qc=rxn.qc.qc))
        template = open(template_file,'r').read()
        template = template.format(label=rxn.instance_name, 
                                   kwargs=kwargs, 
                                   #atom=list(atom),
                                   atom=list(rxn.species.atom),
                                   geom=list([list(gi) for gi in geom]), 
                                   #dummy=dummy,
                                   ppn=rxn.qc.ppn,
                                   qc_command=command,
                                   working_dir=os.getcwd(),
                                   scan=rxn.scan)
    else:
        template_file = pkg_resources.resource_filename('tpl', 'ase_{qc}_ts_end.tpl.py'.format(qc=rxn.qc.qc))
        template = open(template_file,'r').read()
    
        template = template.format(label=rxn.instance_name, 
                                   kwargs=kwargs, 
                                   #atom=list(atom),
                                   atom=list(rxn.species.atom),
                                   geom=list([list(gi) for gi in geom]), 
                                   #dummy=dummy,
                                   ppn=rxn.qc.ppn,
                                   qc_command=command,
                                   working_dir=os.getcwd())
                                   
    with open('{}.py'.format(rxn.instance_name),'w') as f_out:
        f_out.write(template)
    
    step += rxn.qc.submit_qc(rxn.instance_name, singlejob=0)

    return step
