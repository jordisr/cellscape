'''
Load PDB and project atomic coordinates down to 2d
'''

from Bio.PDB import *
import numpy as np

parser = PDBParser()
structure = parser.get_structure('4NOB', '4NOB.pdb')
chain = structure[0]['A']
coords = np.array([list(atom.get_vector()) for atom in chain.get_atoms()])

# as a test, just take x-y and drop z dimension
# output coordinates to csv
pts = coords[:,:2]
np.savetxt('data/test_coords.csv', np.array(pts).T, delimiter=',')
