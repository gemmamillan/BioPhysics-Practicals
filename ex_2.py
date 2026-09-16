#Generate a list of all atoms for a given residue number
#Parameters: PDB file name, Residue number (Including Chain if applicable)

from Bio.PDB import *
import argparse
import numpy as np

parser = argparse.ArgumentParser(
                                 prog='ProgName', 
                                 description='Description of the program'
                                 )
#structure = parser.get_structure("PDB_file_name", "pdb1fat.ent")

parser.add_argument('PDB_file_name',
                    help='Required PDB structure for the program',
                    )
parser.add_argument('residue_number',
                    type = int,
                    help = 'Required residue number')

parser.add_argument('--chain_number')

args = parser.parse_args()

PDBFILE = args.PDB_file_name
residue_num = args.residue_number

PDBparser = PDBParser()

pdbl = PDBList()
downloaded_path = pdbl.retrieve_pdb_file(PDBFILE, pdir='.', file_format='pdb')
st = PDBparser.get_structure(PDBFILE, downloaded_path)

selected = []

if args.chain_number:
    for at in st.get_atoms():
       residue = at.get_parent()
       chain = residue.get_parent()
       if chain.id == args.chain_number:
            selected.append(at)
    print(f'Chain {args.chain_number} from {PDBFILE} has these atoms: {selected}')

    
else: #no chain id was parsed by the user!
    for at in st.get_atoms():
        residue = at.get_parent()
        # print(residue.id)
        # print(type(residue.id[1]))
        if residue.id[1] == residue_num:
            selected.append(at)

    print(f'Residue {residue_num} from {PDBFILE} has these atoms: {selected}')
# for atom in selected:
#     print(f"{atom.get_parent().get_resname()}, {atom.get_parent().id}, {atom.get_name()}, {atom.get_coord()}")