# 5. Generate a list of backbone connectivity (i.e. which residues are linked by ordinary
# peptide bonds).
# Parameters: PDB file name. Optional: Cut-off distance for peptide bonds (defaults to
# # 2.5)
# 5. Ordinary peptide bonds are made between atom C of one residue and atom N of the
# following, Usual distance should be below 2 Å. Following the same approach in exercises 1, or
# 3, find pairs of C-N atoms from different residues that are closer than 2 Å.

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
                    help='Required distance to determine all possible natural peptide bonds')


# Read command line into args
args = parser.parse_args()

cut_distance = 2.5

if args.cut_off_distance:
    set_distance = float(args.cut_off_distance)

    if set_distance > 2.5:
        print('Cut-off distance must be less than 2.5 Å')
        sys.exit(1)
        
    else:
        cut_distance = set_distance

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

#print the variables once assigned
PDBFILE = args.PDB_file_name


# load structure from PDB file
pdbl = PDBList()
downloaded_path = pdbl.retrieve_pdb_file(PDBFILE, pdir='.', file_format='pdb')
st = PDBparser.get_structure(PDBFILE, downloaded_path)

C_atoms = []
N_atoms = []

#Select only CA atoms
for at in st.get_atoms():
    if at.id == 'C':
        C_atoms.append(at)

    elif at.id == 'N':
        N_atoms.append(at)


# Preparing search
nbsearch_n = NeighborSearch(N_atoms)

c_n_pairs = {}

for c_at in C_atoms:
    near_n_atoms = nbsearch_n.search(c_at.get_coord(), cut_distance)

    for n_at in near_n_atoms:
        distance = c_at - n_at
        resc = c_at.get_parent()
        resn = n_at.get_parent()
        pair = (
                f'Carbon residue {resc.get_resname()}{resc.get_id()[1]}: {c_at.get_name()}',
                f'Nitrogen residue {resn.get_resname()}{resn.get_id()[1]}: {n_at.get_name()}'
                )
        c_n_pairs[pair] = round(float(distance), 2)


print(c_n_pairs)
print('Output -> pairs of C and N atoms from different residues closer than the cut-off: { (C atom, N atom) : distance in Å }')