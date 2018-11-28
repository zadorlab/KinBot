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
        'ase_gauss_freq_well.py.tpl',
        'ase_gauss_hir.py.tpl',
        'ase_gauss_irc.py.tpl',
        'ase_gauss_opt_well.py.tpl',
        'ase_gauss_ts_end.py.tpl',
        'ase_gauss_ts_search.py.tpl',
        'ase_nwchem_freq_well.py.tpl',
        'ase_nwchem_irc.py.tpl',
        'ase_nwchem_opt_well.py.tpl',
        'ase_nwchem_ts_end.py.tpl',
        'ase_nwchem_ts_search.py.tpl',
        'ase_nwchem_ts_search_ase_constraints.py.tpl',
        'gauss_freq.tpl',
        'gauss_HO2_Elimination_from_PeroxyRadical_0.tpl',
        'gauss_HO2_Elimination_from_PeroxyRadical_1.tpl',
        'gauss_irc.tpl',
        'gauss_opt_ts.tpl',
        'gauss_opt_well.tpl',
        'gauss_singlept.tpl',
        'gauss_ts_scan_mid.tpl',
        'gauss_ts_scan_start.tpl',
        'gauss_ts_search_end.tpl',
        'gauss_ts_search_mid.tpl',
        'gauss_ts_search_start.tpl',
        'mess_atom.tpl',
        'mess_bimol.tpl',
        'mess_dummy.tpl',
        'mess_fragment.tpl',
        'mess_header.tpl',
        'mess_hinderedrotor.tpl',
        'mess_ts.tpl',
        'mess_tunneling.tpl',
        'mess_well.tpl',
        'molpro.tpl',
        'nwchem.tpl',
        'nwchem_irc.tpl',
        'nwchem_ts_search_end.tpl',
        'nwchem_ts_search_mid.tpl',
        'nwchem_ts_search_start.tpl',
        'pbs.tpl',
        'pbs_mesmer.tpl',
        'pbs_mess.tpl',
        'pbs_python.tpl',
        'slurm_python.tpl']},
    include_package_data=True,
    entry_points={'console_scripts':[
        'kinbot = kinbot.kb:main',
        'pes = kinbot.pes:main',
        ]},
    install_requires=['ase','numpy'],
    
    author="Judit Zador and Ruben Van de Vijver",
    author_email = "jzador@sandia.gov",
    description = "Automatic Potential Energy Surface searches to identify chemical reactions.",
    license = "BSD 3-clause",
    url = "https://github.com/zadorlab/KinBot",
)
