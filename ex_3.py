# Determine all possible hydrogen bonds (Polar atoms at less than 3.5 Å).
# Parameters: PDB file name. Optional: cut-off distance (defaults to 3.5)

#!/usr/bin/env python
#

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
                    help='Required distance to determine all possible hydrogen bonds between polar atom pairs')


# Read command line into args
args = parser.parse_args()

cut_distance = 3.5

if args.cut_off_distance:#(Polar atoms at less than 3.5 Å)
    set_distance = float(args.cut_off_distance)

    if set_distance > 3.5:
        print('Cut-off distance must be less than 3.5 Å')
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

polar_atoms = []

#Select only polar atoms
for at in st.get_atoms():
    if at.id == 'O' or at.id == 'N' or at.id == 'S':#only polar atoms
        polar_atoms.append(at)

# Preparing search
nbsearch = NeighborSearch(polar_atoms)
hydrogen_bonded_atoms = {}

for at1, at2 in nbsearch.search_all(cut_distance):
    d = at1 - at2
    res1 = at1.get_parent()
    res2 = at2.get_parent()
    pair = (
        f'{res1.get_resname()}{res1.get_id()[1]}:{at1.get_name()}',
        f'{res2.get_resname()}{res2.get_id()[1]}:{at2.get_name()}'
    )
    hydrogen_bonded_atoms[pair] = round(float(d), 2)


print(hydrogen_bonded_atoms)
print('Output -> pair of polar atoms that form a hydrogen bond: { (atom1 , atom2) : distance between atoms - Å }')