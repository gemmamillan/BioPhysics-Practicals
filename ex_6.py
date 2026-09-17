# #Id 4, but for disulphide bonds
# 6. Disulphide bonds are formed between S atoms of Cys Residues when they are at the
# appropriate distance (around 1.9 Å). Using the same approach as 1, 3, or 5, find S-S contacts.
# Allow some more distance to access structure variability

from Bio.PDB import *
import argparse
import numpy as np
import sys


PDBparser = PDBParser(PERMISSIVE=1)
parser = argparse.ArgumentParser(
                                 prog='ProgName', 
                                 description='Description of the program'
                                 )

parser.add_argument('PDB_file_name',
                    help='Required PDB structure for the program',
                    )

parser.add_argument('--cut_off_distance',
                    help='Required distance to determine all possible natural di-sulphide bonds')


# Read command line into args
args = parser.parse_args()

cut_distance = 2.05#default proper distance where di-sulphide bonds form

if args.cut_off_distance:
    set_distance = float(args.cut_off_distance)

    if set_distance > 2.5:#Allow some more distance to access structure variability
        print('Cut-off distance must be less than 2.5 Å')
        sys.exit(1)
        
    else:
        cut_distance = set_distance

print('User CLI:')
for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

if args.cut_off_distance == None:
    print('As cut-off distance was set as None,\n default value will be used:', cut_distance)

#print the variables once assigned
PDBFILE = args.PDB_file_name

# load structure from PDB file
pdbl = PDBList()
downloaded_path = pdbl.retrieve_pdb_file(PDBFILE, pdir='.', file_format='pdb')
st = PDBparser.get_structure(PDBFILE, downloaded_path)

S_atoms = []

#Select only CA atoms
for at in st.get_atoms():
    if at.id == 'SG':#SG is the specific sulfur atom id for CYS residues
        residue = at.get_parent()
        if residue.get_resname() == 'CYS':
            S_atoms.append(at)

# Preparing search
if S_atoms == []:
    print('No S atoms were found in this molecule')
    sys.exit(1)

nbsearch = NeighborSearch(S_atoms)

sulfide_bonded_atoms = {}

for at1, at2 in nbsearch.search_all(cut_distance):
    d = at1 - at2
    res1 = at1.get_parent()
    res2 = at2.get_parent()
    pair = (
        f'{res1.get_resname()}{res1.get_id()[1]}: {at1.get_name()}',
        f'{res2.get_resname()}{res2.get_id()[1]}: {at2.get_name()}'
    )
    sulfide_bonded_atoms[pair] = round(float(d), 2)

print(sulfide_bonded_atoms)
print('Output -> pair of sulfur atoms that form a disulfide bond: { (residue1 id1 : sulfur atom forming the disulfide bond, residue2 id2 : sulfur atom) : distance of the bond  - Å }')

