# Print distances between all atom pairs of two given residues
# Parameters: PDB file name, Residue 1, Residue 2
# . Select the desired residues and follow ex_distances_2res.py example

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

parser.add_argument('residue1',
                    help='Required residue to determine all atom distances')

parser.add_argument('residue2',
                    help='Required residue to determine all atom distances')

# Read command line into args
args = parser.parse_args()


print('User CLI:')
for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

#print the variables once assigned
PDBFILE = args.PDB_file_name
residue1 = args.residue1
residue2 = args.residue2

# load structure from PDB file
pdbl = PDBList()
downloaded_path = pdbl.retrieve_pdb_file(PDBFILE, pdir='.', file_format='pdb')
st = PDBparser.get_structure(PDBFILE, downloaded_path)

residues_1 = []
residues_2 = []
#finding the user selected pair of residues
for res in st.get_residues():
    if res.get_resname() == residue1:
        residues_1.append(res)

    elif res.get_resname() == residue2:
        residues_2.append(res)

# #Selection of Residue 10 of Chain A
# res10 = st[0]["A"][10]
# #Selection of Residue 20 of Chain A
# res20 = st[0]["A"][20]

result = {}

for r1 in residues_1:
    for at1 in r1.get_atoms():
        for r2 in residues_2:
            for at2 in r2.get_atoms():
                dist = at1 - at2
                pair = (
                f'{r1.get_resname()}{r1.get_id()[1]}: {at1.get_name()}',
                f'{r2.get_resname()}{r2.get_id()[1]}: {at2.get_name()}'
                )
                result[pair] = round(float(dist), 2)

print(result)
print('Output -> distance between all of the atoms of 2 given residues { (residue1+id: atom+id , residue2+id:atom+id) : distance between atoms - Å }')

