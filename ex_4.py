# Generate a list of all CA atoms of given residue type with coordinates
# Parameters: PDB file name, residue type
# . Iterate over all residues and select those matching the required residue type. For residues
# selected, print the required information

#1. Determine the list of pairs of residues whose CA atoms are closer than a given distance
#Parameters: PDB file name, distance

#!/usr/bin/env python
#
""" Simple program to search contacts """

from Bio.PDB import *
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
import argparse
import numpy as np
import sys


PDBparser = PDBParser(PERMISSIVE=1)
parser = argparse.ArgumentParser(
                                 prog='ProgName', 
                                 description='Description of the program'
                                 )
#structure = parser.get_structure("PDB_file_name", "pdb1fat.ent")

parser.add_argument('PDB_file_name',
                    help='Required PDB structure for the program',
                    )


parser.add_argument('residue_type',
                    help='Required residue specification')


# Read command line into args
args = parser.parse_args()

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

#print the variables once assigned
PDBFILE = args.PDB_file_name
RESIDUE = args.residue_type

# load structure from PDB file
pdbl = PDBList()
downloaded_path = pdbl.retrieve_pdb_file(PDBFILE, pdir='.', file_format='pdb')
st = PDBparser.get_structure(PDBFILE, downloaded_path)

CA_atoms = []
CA_residues = {}

#Select only CA atoms
for at in st.get_atoms():
    if at.id == 'CA':
        CA_atoms.append(at)

if CA_atoms == []:
    print('No CA atoms were found in this molecule')
    sys.exit(1)

for atom in CA_atoms:
    residue = atom.get_parent()
    if residue.get_resname() == RESIDUE:
        store = (f'{residue.get_resname()} {residue.get_id()[1]}, {atom.get_name()}')
        CA_residues[store] = atom.get_coord()

print(CA_residues)
print('Output -> selected residue and its carbon atoms location: { (selected_residue , CA) : [CA atom coordinate array] }')