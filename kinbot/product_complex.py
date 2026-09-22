"""Check optional encounter wells against the final represented products."""
import logging

from kinbot import constants
from kinbot.conformer_counting import writer_members


def represented_ground(species, par, *, fragment=False):
    """Use the same lowest conformer reference as the MESS well/fragment writer."""
    if (par.get('multi_conf_tst') and species.natom > 1
            and not (fragment and species.chemid == 170170000000000000002)):
        # MESS's existing special OH fragment always uses its selected parent.
        records = writer_members(species, par.get('optical_population', 'specified'))
        if records:
            return min(record.zero_energy_hartree for record in records.values())
    return species.energy + species.zpe


def reassess_product_complex(reaction, par):
    """Keep an optional complex only if bound and within the well-mode tolerance."""
    if not reaction.do_vdW:
        return
    complex_species = reaction.irc_prod_opt.species
    fragments = [opt.species for opt in reaction.prod_opt]
    separated = sum(represented_ground(p, par, fragment=True) for p in fragments)
    bound = represented_ground(complex_species, par)
    depth = (separated - bound) * constants.AUtoKCAL
    # L1-only complexes bypass Optimize.compare_structures. Apply the same
    # large-imaginary-mode rejection as for ordinary final product wells.
    lowest_freq = min(complex_species.freq, default=0.)
    imagfreq_threshold = par.get('imagfreq_threshold', 50.)
    frequency_ok = lowest_freq >= -imagfreq_threshold
    reaction.vdW_depth = depth
    reaction.final_vdW_assessment = dict(
        complex_zero_energy_hartree=bound, products_zero_energy_hartree=separated,
        depth_kcal_mol=depth, threshold_kcal_mol=par['vdW_detection'],
        lowest_frequency_cm1=lowest_freq, imagfreq_threshold_cm1=imagfreq_threshold,
        retained=bool(depth > par['vdW_detection'] and frequency_ok))
    if not frequency_ok:
        reaction.do_vdW = False
        logging.getLogger('KinBot').warning(
            '%s: final product-complex frequency %.4f cm-1 is below the allowed '
            '-%.4f cm-1; omitting the optional complex and retaining the TS '
            'and separated products.', reaction.instance_name, lowest_freq, imagfreq_threshold)
    elif depth <= par['vdW_detection']:
        reaction.do_vdW = False
        logging.getLogger('KinBot').info(
            '%s: final product-complex depth %.6f kcal/mol does not exceed %.6f; '
            'omitting the optional complex and retaining the TS and separated products.',
            reaction.instance_name, depth, par['vdW_detection'])
