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

"""
This file is used to install KinBot.

Type 
python setup.py build
python setup.py install
"""

from setuptools import setup, find_packages

setup(
    name = "KinBot",
    version = "2.0",
    packages = find_packages(),
    package_data={'tpl':[
        'arrow.png',
        'ase_gauss_freq_well.tpl.py',
        'ase_gauss_hir.tpl.py',
        'ase_gauss_irc.tpl.py',
        'ase_gauss_opt_well.tpl.py',
        'ase_gauss_ts_end.tpl.py',
        'ase_gauss_ts_search.tpl.py',
        'ase_nwchem_freq_well.tpl.py',
        'ase_nwchem_irc.tpl.py',
        'ase_nwchem_opt_well.tpl.py',
        'ase_gauss_ring_conf.tpl.py',
        'ase_nwchem_ts_end.tpl.py',
        'ase_nwchem_ts_search.tpl.py',
        'ase_nwchem_ts_search_ase_constraints.tpl.py',
        'mess_2tst.tpl',
        'mess_atom.tpl',
        'mess_barrier.tpl',
        'mess_barrierless.tpl',
        'mess_bimol.tpl',
        'mess_core_rr.tpl',
        'mess_dummy.tpl',
        'mess_fragment.tpl',
        'mess_fragment_OH.tpl',
        'mess_pstfragment.tpl',
        'mess_header.tpl',
        'mess_hinderedrotor.tpl',
        'mess_outerrrho.tpl',
        'mess_pst.tpl',
        'mess_rrho.tpl',
        'mess_termol.tpl',
        'mess_tunneling.tpl',
        'mess_variational.tpl',
        'mess_well.tpl',
        'molpro.tpl',
        'pbs.tpl',
        'pbs_molpro.tpl',
        'pbs_mesmer.tpl',
        'pbs_mess.tpl',
        'pbs_python.tpl',
        'pesviewer.inp.tpl',
        'slurm.tpl',
        'slurm_molpro.tpl',
        'slurm_mess.tpl',
        'slurm_mesmer.tpl',
        'slurm_python.tpl']},
    include_package_data=True,
    entry_points={'console_scripts':[
        'kinbot = kinbot.kb:main',
        'pes = kinbot.pes:main',
        ]},
    install_requires=['ase','numpy','networkx'],
    
    author="Judit Zador and Ruben Van de Vijver",
    author_email = "jzador@sandia.gov",
    description = "Automatic Potential Energy Surface searches to identify chemical reactions.",
    license = "BSD 3-clause",
    url = "https://github.com/zadorlab/KinBot",
)
