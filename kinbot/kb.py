import sys
import logging
import datetime
import copy

from kinbot import license_message
from kinbot import postprocess
from kinbot.parameters import Parameters
from kinbot.mess import MESS
from kinbot.optimize import Optimize
from kinbot.reaction_finder import ReactionFinder
from kinbot.reaction_finder_bimol import ReactionFinderBimol
from kinbot.reaction_generator import ReactionGenerator
from kinbot.stationary_pt import StationaryPoint
from kinbot.qc import QuantumChemistry
from kinbot.utils import make_dirs, clean_files
from kinbot.config_log import config_log


def main():
    if sys.version_info.major < 3:
        print(f'KinBot only runs with python 3.8 or higher. You have python {sys.version_info.major}.{sys.version_info.minor}. Bye!')
        sys.exit(-1)
    elif sys.version_info.minor < 8:
        print(f'KinBot only runs with python 3.8 or higher. You have python {sys.version_info.major}.{sys.version_info.minor}. Bye!')
        sys.exit(-1)

    try:
        input_file = sys.argv[1]
    except IndexError:
        print('To use KinBot, supply one argument being the input file!')
        sys.exit(-1)

    print(license_message.message)

    # initialize the parameters for this run
    masterpar = Parameters(input_file)
    par = masterpar.par
    input_file = masterpar.input_file
    # set up the logging environment
    if par['verbose']:
        logger = config_log('KinBot', 'debug')
    else:
        logger = config_log('KinBot')

    # write the license message and the parameters to the log file
    logger.info(license_message.message)
    logger.info('Input parameters')
    par_str = "\n\t".join([str(p) + ": " + str(par[p])for p in par])

    logger.info(par_str)
    # time stamp of the KinBot start
    logger.info('Starting KinBot')

    make_dirs(par)

    if par['bimol'] == 0:
        # initialize the reactant
        well0 = StationaryPoint('well0',
                                par['charge'],
                                par['mult'],
                                smiles=par['smiles'],
                                structure=par['structure'])
        well0.short_name = 'w1'
        # write the initial reactant geometry to a file for visualization
        with open('geometry.xyz', 'w') as geom_out:
            geom_out.write('{}\n\n'.format(well0.natom))
            for i, at in enumerate(well0.atom):
                x, y, z = well0.geom[i]
                geom_out.write('{} {:.6f} {:.6f} {:.6f}\n'.format(at, x, y, z))
            geom_out.write('\n\n')

        # characterize the initial reactant
        well0.characterize()
        well0.name = str(well0.chemid)
        start_name = well0.name
        if well0.name in par['skip_chemids']:
            logger.info('This chemid is skipped, nothing to do here')
            logger.info('Finished KinBot at {}'.format(datetime.datetime.now()))
            print("Done!")
            return

        # initialize the qc instance
        qc = QuantumChemistry(par)

        clean_files()

        # start the initial optimization of the reactant
        logger.info('Starting optimization of initial well...')
        qc.qc_opt(well0, well0.geom)
        err, well0.geom = qc.get_qc_geom(str(well0.chemid) + '_well',
                                         well0.natom, wait=1)
        logger.debug(f'Initial well opt error {err}.')
        err, well0.freq = qc.get_qc_freq(str(well0.chemid) + '_well',
                                         well0.natom, wait=1)
        logger.debug(f'Initial well freq error {err}, frequencies are {well0.freq}.')
        if err < 0:
            logger.error('Error with initial structure optimization.')
            return
        if well0.freq[0] <= 0:
            logger.warning(f'First frequency is {well0.freq[0]} for initial structure.')
            well0.freq[0] *= -1.
        if well0.freq[1] <= 0:
            logger.error(f'Second frequency is {well0.freq[1]} for initial structure.')
            return

        # characterize again and look for differences
        well0 = StationaryPoint('well0',
                                par['charge'],
                                par['mult'],
                                atom=copy.deepcopy(well0.atom),
                                geom=copy.deepcopy(well0.geom))
        well0.short_name = 'w1'
        well0.characterize()
        well0.name = str(well0.chemid)
        if well0.name != start_name:
            logger.error('The first well optimized to a structure different from the input.')
            return

        # do an MP2 optimization of the reactant,
        # to compare some scan barrier heigths to
        if par['families'] == ['all'] or \
                'birad_recombination_R' in par['families'] or \
                'r12_cycloaddition' in par['families'] or \
                'r14_birad_scission' in par['families'] or \
                'R_Addition_MultipleBond' in par['families'] or \
                (par['skip_families'] != ['none'] and
                ('birad_recombination_R' not in par['skip_families'] or
                'r12_cycloaddition' not in par['skip_families'] or
                'r14_birad_scission' not in par['skip_families'] or
                'R_Addition_MultipleBond' not in par['skip_families'])) or \
                par['reaction_search'] == 0:
            logger.debug('Starting MP2 optimization of initial well...')
            qc.qc_opt(well0, well0.geom, mp2=1)
            err, geom = qc.get_qc_geom(str(well0.chemid) + '_well_mp2', well0.natom, 1)

        # comparison for barrierless scan
        if par['barrierless_saddle']:
            logger.debug('Optimization of intial well for barrierless at {}/{}'.
                    format(par['barrierless_saddle_method'], par['barrierless_saddle_basis']))
            qc.qc_opt(well0, well0.geom, bls=1)
            err, geom = qc.get_qc_geom(str(well0.chemid) + '_well_bls', well0.natom, 1)

        # characterize again and look for differences
        well0.characterize()
        well0.name = str(well0.chemid)

        err, well0.energy = qc.get_qc_energy(str(well0.chemid) + '_well', 1)
        err, well0.zpe = qc.get_qc_zpe(str(well0.chemid) + '_well', 1)
        # to save the starting energy to make thresholds consistent
        well0.start_energy = well0.energy
        well0.start_zpe = well0.zpe

        well_opt = Optimize(well0, par, qc, wait=1)
        well_opt.do_optimization()
        if well_opt.shigh == -999:
            logger.error('Error with high level optimization of initial structure.')
            return

        # if par['pes']:
        #    filecopying.copy_to_database_folder(well0.chemid, well0.chemid, qc)

        if par['reaction_search'] == 1:
            logger.info('\tStarting reaction search...')
            rf = ReactionFinder(well0, par, qc)
            rf.find_reactions()
            rg = ReactionGenerator(well0, par, qc, input_file)
            rg.generate()

    # BIMOLECULAR REACTANTS
    elif par['bimol'] == 1:
        fragments = {'frag_a': None, 'frag_b': None} 
        # initialize the reacting fragments
        charge = 0
        for ii, frag in enumerate(fragments):
            fragments[frag] = StationaryPoint(frag,
                                              par['charge'][ii],
                                              par['mult'][ii],
                                              structure=par['structure'][ii])
            fragments[frag].characterize()
            fragments[frag].name = str(fragments[frag].chemid)
            charge += par['charge'][ii]
         
        # this is just formally a well
        well0 = StationaryPoint('bimolecular',
                                par['charge'],
                                par['mult'][2],
                                smiles=par['smiles'],
                                structure=par['structure'],
                                fragA=fragments['frag_a'],
                                fragB=fragments['frag_b'],
                                )
        well0.characterize()
        well0.short_name = 'w1'
        well0.chemid = '_'.join(sorted([fragments['frag_a'].name, fragments['frag_b'].name]))
        well0.name = str(well0.chemid)
 
        qc = QuantumChemistry(par)

        logger.info('\tStarting optimization of fragments...')
        for frag in fragments.values():
            qc.qc_opt(frag, frag.geom)
            err, frag.geom = qc.get_qc_geom(str(frag.chemid) + '_well',
                                            frag.natom, wait=1)
            if err < 0:
                logger.error(f'Error with initial structure optimization of {frag.name}.')
                return
            err, frag.freq = qc.get_qc_freq(str(frag.chemid) + '_well',
                                            frag.natom, wait=1)
            if frag.freq[0] <= 0. and frag.freq[0] >= -20.:
                logger.warning(f'Found imaginary frequency {frag.freq[0]} for {frag.name}. It is flipped.')
                frag.freq[0] *= -1.
            elif frag.freq[0] < -20.:
                logger.error(f'Found imaginary frequency {frag.freq[0]} for {frag.name}.')
                return
            if frag.freq[1] <= 0:
                logger.error(f'Found two imaginary frequencies {frag.freq[1]} for {frag.name}.')
                return
            # characterize again and look for differences
            frag.characterize()
            if frag.name != str(frag.chemid):
                logger.error(f'Reactant {frag} optimized to a structure different from the input.')
                return

        well0.energy = 0.
        well0.zpe = 0.
        well0.fragA = fragments['frag_a']  # update
        well0.fragB = fragments['frag_b']  # update
        well0.structure = [well0.fragA.structure, well0.fragB.structure]  # update, used in reactions
        for frag in fragments.values():
            err, frag.energy = qc.get_qc_energy(str(frag.chemid) + '_well', 1)
            err, frag.zpe = qc.get_qc_zpe(str(frag.chemid) + '_well', 1)
            well0.energy += frag.energy
            well0.zpe += frag.zpe

            frag_opt = Optimize(frag, par, qc, wait=1)
            frag_opt.do_optimization()
            if frag_opt.shigh == -999:
                logger.error(f'Error with high level optimization for {frag.name}.')
                return

        well0.fragA = fragments['frag_a']  # update
        well0.fragB = fragments['frag_b']  # update
        well0.energy = 0.
        well0.zpe = 0.
        for frag in fragments.values():
            well0.energy += frag.energy
            well0.zpe += frag.zpe
        # if par['pes']:
        #    filecopying.copy_to_database_folder(well0.chemid, well0.chemid, qc)

        if par['reaction_search'] == 1:
            logger.info('\tStarting bimolecular reaction search...')
            rf = ReactionFinderBimol(well0, par, qc)
            rf.find_reactions()
            rg = ReactionGenerator(well0, par, qc, input_file)
            rg.generate()

    if par['me'] > 0:  # it will be 2 for kinbots when the mess file is needed but not run
        mess = MESS(par, well0)
        mess.write_input(qc)

        if par['me'] == 1:
            logger.info('\tStarting Master Equation calculations')
            if par['me_code'] == 'mess':
                mess.run()

    postprocess.create_summary_file(well0, qc, par)
    postprocess.createPESViewerInput(well0, qc, par)
    postprocess.creatMLInput(well0, qc, par)

    logger.info('Finished KinBot at {}'.format(datetime.datetime.now()))
    try:
        print("Done!")
    except OSError:
        pass
    return


if __name__ == "__main__":
    main()
