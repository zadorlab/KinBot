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
