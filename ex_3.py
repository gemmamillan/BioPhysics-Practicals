# Determine all possible hydrogen bonds (Polar atoms at less than 3.5 Å,
# filtered to pairs where one atom can actually act as a donor (has an H
# to give) and the other as an acceptor (has a lone pair to receive it)).
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



BACKBONE_NO_DONOR_RESIDUES = {'PRO'}

SIDECHAIN_ROLES = {
    ('ARG', 'NE'):  {'donor'},
    ('ARG', 'NH1'): {'donor'},
    ('ARG', 'NH2'): {'donor'},
    ('ASN', 'ND2'): {'donor'},
    ('ASN', 'OD1'): {'acceptor'},
    ('ASP', 'OD1'): {'acceptor'},
    ('ASP', 'OD2'): {'acceptor'},
    ('CYS', 'SG'):  {'donor', 'acceptor'},   # thiol, weak H-bonder
    ('GLN', 'NE2'): {'donor'},
    ('GLN', 'OE1'): {'acceptor'},
    ('GLU', 'OE1'): {'acceptor'},
    ('GLU', 'OE2'): {'acceptor'},
    ('HIS', 'ND1'): {'donor', 'acceptor'},   # depends on tautomer/protonation
    ('HIS', 'NE2'): {'donor', 'acceptor'},
    ('LYS', 'NZ'):  {'donor'},               # protonated amine
    ('MET', 'SD'):  {'acceptor'},            # weak acceptor, no H
    ('SER', 'OG'):  {'donor', 'acceptor'},   # hydroxyl
    ('THR', 'OG1'): {'donor', 'acceptor'},   # hydroxyl
    ('TRP', 'NE1'): {'donor'},
    ('TYR', 'OH'):  {'donor', 'acceptor'},   # hydroxyl
}


def get_atom_role(atom):
    """Return the set of roles ({'donor'}, {'acceptor'}, both, or
    empty) this atom can play in a hydrogen bond, based on the
    residue it belongs to and its atom name."""
    residue = atom.get_parent()
    resname = residue.get_resname()
    name = atom.get_name()

    # Water: the O can both donate (via either H) and accept.
    if resname == 'HOH':
        return {'donor', 'acceptor'}

    # Backbone amide nitrogen: donor, except proline (tertiary N, no H).
    if name == 'N':
        if resname in BACKBONE_NO_DONOR_RESIDUES:
            return set()
        return {'donor'}

    # Backbone carbonyl oxygen (and terminal carboxylate O): acceptor only.
    if name in ('O', 'OXT'):
        return {'acceptor'}

    # Side-chain lookup table.
    if (resname, name) in SIDECHAIN_ROLES:
        return SIDECHAIN_ROLES[(resname, name)]

    # Fallback by element, for anything not explicitly covered.
    element = (atom.element or name[0]).strip()
    if element == 'N':
        return {'donor'}
    if element == 'O':
        return {'acceptor'}
    if element == 'S':
        return {'acceptor'}
    return set()


def is_real_hbond(role1, role2):
    """True only if one side can donate and the other can accept."""
    if not role1 or not role2:
        return False
    return ('donor' in role1 and 'acceptor' in role2) or \
           ('acceptor' in role1 and 'donor' in role2)

# Read command line into args
args = parser.parse_args()

cut_distance = 3.5

if args.cut_off_distance:  # (Polar atoms at less than 3.5 Å)
    set_distance = float(args.cut_off_distance)

    if set_distance > 3.5:
        print('Cut-off distance must be less than 3.5 Å')
        sys.exit(1)

    else:
        cut_distance = set_distance

for k, v in vars(args).items():
    print('{:10}:'.format(k), v)

# print the variables once assigned
PDBFILE = args.PDB_file_name

# load structure from PDB file
pdbl = PDBList()
downloaded_path = pdbl.retrieve_pdb_file(PDBFILE, pdir='.', file_format='pdb')
st = PDBparser.get_structure(PDBFILE, downloaded_path)

polar_atoms = []

# Select only polar atoms
for at in st.get_atoms():
    if at.id == 'O' or at.id == 'N' or at.id == 'S':  # only polar atoms
        polar_atoms.append(at)

# Preparing search
nbsearch = NeighborSearch(polar_atoms)
hydrogen_bonded_atoms = {}

for at1, at2 in nbsearch.search_all(cut_distance):

    role1 = get_atom_role(at1)
    role2 = get_atom_role(at2)

    # Keep the pair only if a donor-acceptor relationship is chemically possible.
    if not is_real_hbond(role1, role2):
        continue

    d = at1 - at2
    res1 = at1.get_parent()
    res2 = at2.get_parent()
    pair = (
        f'{res1.get_resname()}{res1.get_id()[1]}: {at1.get_name()}',
        f'{res2.get_resname()}{res2.get_id()[1]}: {at2.get_name()}'
    )
    hydrogen_bonded_atoms[pair] = round(float(d), 2)


print(hydrogen_bonded_atoms)
print('Output -> pair of polar atoms that form a hydrogen bond: { (atom1 , atom2) : distance between atoms - Å }')