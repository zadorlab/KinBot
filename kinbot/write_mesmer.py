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
import os,sys
import xml.etree.cElementTree as ET
import xml.dom.minidom as minidom
import random

def write_mesmer_input(species,barriers,products):
    root = ET.Element(  'me:mesmer',{'xmlns':'http://www.xml-cml.org/schema',
                        'xmlns:me':'http://www.chem.leeds.ac.uk/mesmer',
                        'xmlns:xsi':'http://www.w3.org/2001/XMLSchema-instance'})
    
    title = ET.SubElement(root,'me:title').text = 'species.chemid'
    mollist = ET.SubElement(root,'moleculeList')
    
    #write the initial species
    atom = ['C','C','C','H','H','H','H','H','H']
    natom = len(atom)
    rad = [0 for ai in atom]
    charge = 0
    addMolecule(mollist, species, atom, natom, rad, charge)
    
    #Todo: write the products and tss to the mollist
    
    reaclist = ET.SubElement(root,'reactionList')
    #write the reactions
    for index, instance in enumerate(species.reac_inst):
        addReaction(reaclist, species, index, instance)
    
    st = ET.tostring(root,'utf-8')
    st = minidom.parseString(st)
    fout = open('test.xml','w')
    fout.write(st.toprettyxml(indent = ' '))
    fout.close()
    #write st.toprettyxml(indent = ' ')
    #tree.write('test.xml', encoding='utf-8',xml_declaration=True)

def addReaction(reaclist, species, index, instance):
    a = 1

def addMolecule(mollist,mol, atom, natom, rad, charge):
    geom = []
    for i,at in enumerate(atom):
        geom.append([random.uniform(-3.,3.), random.uniform(-3.,3.), random.uniform(-3.,3.)])
    bond = [
    [0,2,0,1,1,0,0,0,0],
    [2,0,1,0,0,1,0,0,0],
    [0,1,0,0,0,0,1,1,1],
    [1,0,0,0,0,0,0,0,0],
    [1,0,0,0,0,0,0,0,0],
    [0,1,0,0,0,0,0,0,0],
    [0,0,1,0,0,0,0,0,0],
    [0,0,1,0,0,0,0,0,0],
    [0,0,1,0,0,0,0,0,0]
    ]
    molecule = ET.SubElement(mollist, 'molecule', {'id':'species.chemid','spinMultiplicity':'{}'.format(sum(rad))})
    atomarray = ET.SubElement(molecule, 'atomArray')
    for i,at in enumerate(atom):
        args = {'id':'a{}'.format(i+1)}
        args['elementType'] = at
        args['x3'] = '{:.8f}'.format(geom[i][0])
        args['y3'] = '{:.8f}'.format(geom[i][1])
        args['z3'] = '{:.8f}'.format(geom[i][2])
        at = ET.SubElement(atomarray, 'atom', args)
    bond_id = 1
    bondarray = ET.SubElement(molecule, 'bondArray')
    for i in range(len(atom)-1):
        for j in range(i+1,len(atom)):
            if bond[i][j] > 0:
                args = {'id':'b{}'.format(bond_id)}
                args['atomRefs2']="a{} a{}".format(i+1,j+1)
                args['order']="{}".format(bond[i][j])
                b = ET.SubElement(bondarray,'bond',args)
                bond_id += 1
    propertylist = ET.SubElement(molecule, 'propertyList')
    
    #add the zpe
    property = ET.SubElement(propertylist, 'property', {'dictRef':'me:ZPE'})
    scalar = ET.SubElement(property, 'scalar', {'units':'cm-1'}).text = str(15.5)
    

def indent(elem, level=0):
    i = "\n" + level*"  "
    j = "\n" + (level-1)*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for subelem in elem:
            indent(subelem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = j
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = j
    return elem

species = 'a'
barriers = ['b1','b2']
products = [['p1','p2'],['p3']]
write_mesmer_input(species,barriers,products)