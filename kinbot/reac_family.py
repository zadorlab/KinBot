import os
import numpy as np

from kinbot import kb_path
from kinbot import modify_geom
from kinbot import geometry
from kinbot.reactions.reac_abstraction import abstraction_align


def carry_out_reaction(rxn, step, command, bimol=0):
    """
    Verify what has been done and what needs to be done
    skip: boolean which tells to skip the first 12 steps in case of an instance 
        shorter than 4
    scan: boolean which tells if this is part of an energy scan along a bond 
        length coordinate
    """
    ts = True

    if step > 0:
        status = rxn.qc.check_qc(rxn.instance_name)
        if status != 'normal' and status != 'error':
            return step
  
    kwargs = rxn.qc.get_qc_arguments(rxn.instance_name, rxn.species.mult, 
                                     rxn.species.charge, ts=ts, step=step, 
                                     max_step=rxn.max_step, scan=rxn.scan)
    if step == 0:
        if rxn.qc.is_in_database(rxn.instance_name):
            if rxn.qc.check_qc(rxn.instance_name) == 'normal':  # log file is present and is in the db
                err, freq = rxn.qc.get_qc_freq(rxn.instance_name, rxn.species.natom)
                if err == 0 and len(freq) > 0.:  # only final calculations have frequencies
                    err, geom = rxn.qc.get_qc_geom(rxn.instance_name, rxn.species.natom)
                    step = rxn.max_step + 1  # this shortcuts the search, jumps to the end
                    return step
            if rxn.qc.check_qc(rxn.instance_name) == 'error':  # log file is not present
                err, geom = rxn.qc.get_qc_geom(rxn.instance_name, rxn.species.natom, allow_error=1)
                if np.sum(geom) == 0:
                    return -1  # we don't want this to be repeated

        if rxn.skip and len(rxn.instance) < 4:
            step = 12
        geom = rxn.species.geom
        if bimol:
            if rxn.family_name == 'abstraction':
                # gives the reactant and product geometry guesses
                geom, _, _ = abstraction_align(rxn.species.geom, rxn.instance, rxn.species.atom, rxn.species.fragA.natom)

    elif step == rxn.max_step and rxn.scan:
        err, geom = rxn.qc.get_qc_geom(rxn.instance_name, rxn.species.natom, 
                                       allow_error=1, previous=1)
    else:
        err, geom = rxn.qc.get_qc_geom(rxn.instance_name, rxn.species.natom, allow_error=1)
        if bimol:
            if rxn.family_name == 'abstraction':
                # gives the reactant and product geometry guesses
                _, geom_prod, geom_ts = abstraction_align(geom, rxn.instance, rxn.species.atom, rxn.species.fragA.natom)

    step, fix, change, release = rxn.get_constraints(step, geom)

    if step > rxn.max_step:
        return step

    # apply the geometry changes here and fix the coordinates that changed
    change_starting_zero = []
    for c in change:
        if c[0] == 'L':
            c_new = [ci - 1 for ci in c[1:-1]]
            c_new = ['L'] + c_new
        else:
            c_new = [ci - 1 for ci in c[:-1]]
        c_new.append(c[-1])
        change_starting_zero.append(c_new)

    if len(change_starting_zero) > 0 and "frozen" not in rxn.instance_name:
        success, geom = modify_geom.modify_coordinates(rxn.species, 
                                                       rxn.instance_name, geom, 
                                                       change_starting_zero,
                                                       rxn.species.bond,
                                                       write_files=0)
        for c in change:
            fix.append(c[:-1])
        change = []
    elif "frozen" in rxn.instance_name:
        if step != 0:
            tmp_species = rxn.get_frozen_species(distance=rxn.scan_list[step])
            geom = tmp_species.geom
        else:
            geom = rxn.species.geom
        for c in change:
            fix.append(c[:-1])
        change = []

    if rxn.qc.qc == 'gauss' or (rxn.qc.qc == 'nn_pes' and step < rxn.max_step):
        code = 'gaussian'
        Code = 'Gaussian'
        kwargs['addsec'] = ''
        if (not bimol or step == 0) and not rxn.qc.use_sella:
            # here addsec contains the constraints
            for fixi in fix:
                kwargs['addsec'] += f"{' '.join(str(f) for f in fixi)} F\n"
            for chi in change:
                kwargs['addsec'] += f"{' '.join(str(ch) for ch in chi)} F\n"
            for reli in release:
                kwargs['addsec'] += f"{' '.join(str(rel) for rel in reli)} A\n"
        elif bimol and step == 1 and not rxn.qc.use_sella:
            kwargs['addsec'] = f'{rxn.instance[0] + 1} {rxn.instance[2] + 1}\n\n'
            # here addsec needs to contain the product and ts geometries and all the rest of the fluff
            kwargs['addsec'] += f'product geometry guess\n\n{rxn.species.charge} {rxn.species.mult}\n'
            for ii, at in enumerate(rxn.species.atom):
                kwargs['addsec'] += f'{at} {geom_prod[ii][0]} {geom_prod[ii][1]} {geom_prod[ii][2]}\n'
            kwargs['addsec'] += f'\n{rxn.instance[0] + 1} {rxn.instance[2] + 1}\n\n'
            kwargs['addsec'] += f'ts geometry guess\n\n{rxn.species.charge} {rxn.species.mult}\n'
            for ii, at in enumerate(rxn.species.atom):
                kwargs['addsec'] += f'{at} {geom_ts[ii][0]} {geom_ts[ii][1]} {geom_ts[ii][2]}\n'
            kwargs['addsec'] += f'\n{rxn.instance[0] + 1} {rxn.instance[2] + 1}\n\n'
    elif rxn.qc.qc == 'qchem':
        code = 'qchem'
        Code = 'QChem'
        if (not bimol or step == 0) and step < rxn.max_step:
            kwargs['addsec'] = '$opt\nCONSTRAINT\n'
            for fixi in fix:
                if len(fixi) == 2:
                    fix_type = 'stre'
                    val = np.linalg.norm(geom[fixi[0] - 1] - geom[fixi[1] - 1])
                elif len(fixi) == 3:
                    fix_type = 'bend'
                    val = geometry.calc_angle(geom[fixi[0]-1], geom[fixi[1]-1], geom[fixi[2]-1]) * 180 / np.pi
                elif len(fixi) == 4:
                    fix_type = 'tors'
                    val = geometry.calc_dihedral(geom[fixi[0]-1], geom[fixi[1]-1], geom[fixi[2]-1],
                                                 geom[fixi[3]-1])[0]
                    kwargs['addsec'] += f"{fix_type} {' '.join(str(f) for f in fixi)} {val}\n"
            for chi in change:
                dist = np.linalg.norm(geom[chi[0] - 1] - geom[chi[1] - 1])
                kwargs['addsec'] += f"{' '.join(str(ch) for ch in chi)} {dist}\n"
            for reli in release:
                dist = np.linalg.norm(geom[reli][0] - geom[reli][1])
                kwargs['addsec'] += f"{' '.join(str(rel) for rel in reli)} {dist}\n"
            kwargs['addsec'] += 'ENDCONSTRAINT\n$end\n'
        elif bimol and step == 1:
            raise NotImplementedError('Bimolecular reactions are not yet '
                                      'implemented in QChem')
    elif rxn.qc.qc == 'nn_pes' and step >= rxn.max_step:
        code = 'nn_pes'
        Code = 'Nn_surr'

    if step < rxn.max_step:
        if rxn.qc.use_sella:
            kwargs.pop('addsec', None)
            kwargs.pop('opt', None)
            template_file = f'{kb_path}/tpl/ase_sella_ts_search.tpl.py'
        else:
            template_file = f'{kb_path}/tpl/ase_{rxn.qc.qc}_ts_search.tpl.py'
        template = open(template_file,'r').read()
        template = template.format(label=rxn.instance_name, 
                                kwargs=kwargs, 
                                atom=list(rxn.species.atom),
                                geom=list([list(gi) for gi in geom]),
                                bimol=bimol,
                                ppn=rxn.qc.ppn,
                                qc_command=command,
                                working_dir=os.getcwd(),
                                scan=rxn.scan,
                                code=code,  # Sella
                                Code=Code,  # Sella
                                fix=fix,  # Sella
                                sella_kwargs=rxn.par['sella_kwargs']  # Sella
                                )
    else:
        if rxn.par['calc_kwargs']:
            kwargs = rxn.qc.merge_kwargs(kwargs)
        if rxn.qc.use_sella:
            kwargs.pop('freq', None)
            kwargs.pop('opt', None)
            template_file = f'{kb_path}/tpl/ase_sella_ts_end.tpl.py'
        else:
            template_file = f'{kb_path}/tpl/ase_{rxn.qc.qc}_ts_end.tpl.py'
        template = open(template_file, 'r').read()
    
        template = template.format(label=rxn.instance_name, 
                                   kwargs=kwargs,
                                   atom=list(rxn.species.atom),
                                   geom=list([list(gi) for gi in geom]),
                                   ppn=rxn.qc.ppn,
                                   qc_command=command,
                                   working_dir=os.getcwd(),
                                   code=code,  # Sella
                                   Code=Code,  # Sella
                                   sella_kwargs=rxn.par['sella_kwargs']  # Sella
                                   )
                                   
    with open('{}.py'.format(rxn.instance_name),'w') as f_out:
        f_out.write(template)

    step += rxn.qc.submit_qc(rxn.instance_name, singlejob=0, 
                             jobtype=kwargs.pop('method', None))

    return step
