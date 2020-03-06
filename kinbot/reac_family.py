import os
import numpy as np
import copy
import time
import pkg_resources

from kinbot import modify_geom

def carry_out_reaction(rxn, step, command):
    """
    Verify what has been done and what needs to be done
    skip: boolean which tells to skip the first 12 steps in case of an instance shorter than 4
    scan: boolean which tells if this is part of an energy scan along a bond length coordinate
    """
    if step > 0:
        status = rxn.qc.check_qc(rxn.instance_name)
        if status != 'normal' and status != 'error': return step

    kwargs = rxn.qc.get_qc_arguments(   rxn.instance_name, rxn.species.mult, rxn.species.charge, ts=1,
                                        step = step, max_step=rxn.max_step, scan = rxn.scan)

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
    else:
        err, geom = rxn.qc.get_qc_geom(rxn.instance_name, rxn.species.natom, allow_error = 1)

    step, fix, change, release = rxn.get_constraints(step, geom)

    if step > rxn.max_step:
        return step

    #apply the geometry changes here and fix the coordinates that changed
    change_starting_zero = []
    for c in change:
        c_new = [ci - 1 for ci in c[:-1]]
        c_new.append(c[-1])
        change_starting_zero.append(c_new)
    if len(change_starting_zero) >0 :
        success, geom = modify_geom.modify_coordinates(rxn.species, rxn.instance_name, geom, change_starting_zero, rxn.species.bond)
        for c in change:
            fix.append(c[:-1])
        change = []

    kwargs['fix'] = fix
    kwargs['change'] = change
    kwargs['release'] = release

    if step < rxn.max_step:
        template_file = pkg_resources.resource_filename('tpl', 'ase_{qc}_ts_search.tpl.py'.format(qc=rxn.qc.qc))
        template = open(template_file,'r').read()
    else:
        template_file = pkg_resources.resource_filename('tpl', 'ase_{qc}_ts_end.tpl.py'.format(qc=rxn.qc.qc))
        template = open(template_file,'r').read()
    
    template = template.format(label=rxn.instance_name, 
                               kwargs=kwargs, 
                               atom=list(rxn.species.atom), 
                               geom=list([list(gi) for gi in geom]), 
                               ppn=rxn.qc.ppn,
                               qc_command=command,
                               working_dir=os.getcwd())

    f_out = open('{}.py'.format(rxn.instance_name),'w')
    f_out.write(template)
    f_out.close()
    
    step += rxn.qc.submit_qc(rxn.instance_name, 0)

    return step
