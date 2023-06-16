> <comments>
Copyright 2018 National Technology & Engineering Solutions of Sandia,
LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, 
the U.S. Government retains certain rights in this software.

This comment is not interpreted, so store any extra info here.
Keywords are case insensitive. Look at the help below.
IMPORTANT: avoid the use of '2d' and '3d' in the names of species, transition states and reactions
(these strings are employed when generating the 2d and 3d files of the molecules)
If you want to use 3D coordinates, store them in a xyz/ directory in the same directory as the python script

> <id> {id}

> <options> 
title              0         # print a title (1) or not (0) 
units              kcal/mol  # energy units
use_xyz            1         # use xyz, put 0  to switch off
rescale            0         # no rescale , put the well or bimolecular name here to rescale to that value
fh                 9.        # figure height
fw                 18.       # figure width
margin             0.2       # margin fraction on the x and y axis
dpi                120       # dpi of the molecule figures
save               0         # does the plot need to be saved (1) or displayed (0)
write_ts_values    1         # booleans tell if the ts energy values should be written
write_well_values  1         # booleans tell if the well and bimolecular energy values should be written
bimol_color        red       # color of the energy values for the bimolecular products
well_color         blue      # color of the energy values of the wells
ts_color           green     # color or the energy values of the ts, put to 'none' to use same color as line
show_images        1         # boolean tells whether the molecule images should be shown on the graph
rdkit4depict       1         # boolean that specifies which code was used for the 2D depiction
axes_size          10        # font size of the axes
text_size          10        # font size of the energy values on the graph
graph_edge_color   black     # color of graph edge, if set to 'energy', will be scaled accordingly

> <wells> 
{wells}

> <bimolec> 
{bimolecs}

> <ts> 
{ts}

> <barrierless> 
{barrierless}

> <help>
File follows the rules of SD file format for keywords. Keywords are case
insensitive when parsed.
Keywords:
units: units of the energies supplied above

usexyz: use the xyz coordinates of all the species and render a 2D/3D depiction

rescale: energies are rescaled relative to the energy of the species given here 

wells: all the wells of the PES, separated by lines
each line contains the name, the energy, and optionally the smiles

bimolec: all the bimolecular products of the PES, separated by lines
each line contains the name, the energy, and optionally the smiles of both bimolecular products

ts: all the transition states of the PES, separated by lines
each line contains the name, the energy, and the names of the reactant and product

barrierless: all the barrierless reactions of the PES, separated by lines
each line contains the name and the names of the reactant and product
