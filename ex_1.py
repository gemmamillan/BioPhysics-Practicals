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


parser.add_argument('distance',
                    help='Required distance to determine which pair of residues CA atoms are closer')


# Read command line into args
args = parser.parse_args()

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

#print the variables once assigned
PDBFILE = args.PDB_file_name
DISTANCE = float(args.distance)


# load structure from PDB file
pdbl = PDBList()
downloaded_path = pdbl.retrieve_pdb_file(PDBFILE, pdir='.', file_format='pdb')
st = PDBparser.get_structure(PDBFILE, downloaded_path)

CA_atoms = []

#Select only CA atoms
for at in st.get_atoms():
    if at.id == 'CA':
        CA_atoms.append(at)

# Preparing search

if CA_atoms == []:
    print('No CA atoms were found in this molecule')
    sys.exit(1)

nbsearch = NeighborSearch(CA_atoms)

print("Performing Neighbour search, with max distance set to:", DISTANCE, '...')

#Searching for contacts under HBLNK
ncontact = 1
residue_pairs = {}

for at1, at2 in nbsearch.search_all(DISTANCE):
    
    residues = (at1.get_parent().get_resname(),at2.get_parent().get_resname())

    distance = at1 - at2

    residue_pairs[residues] = (at1.get_serial_number(),at2.get_serial_number(),round(float(distance),2))
    ncontact += 1

if residue_pairs == {}:
    print('No CA atoms were closer than the set distance:', DISTANCE)
    sys.exit(1)

print(residue_pairs)
print('Output -> pair of residues whose carbon atoms are closer than the set distance: { (residue1 , residue2) : (CA atom serial number residue1, CA atom serial number residue2, distance between CA atoms - Å) }')