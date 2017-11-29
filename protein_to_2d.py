'''
Generate 2D space-filling outline from a PDB structure.
Large gaps in the structure may produce errors. Domains are not segmented.
'''

from Bio.PDB import *
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import shapely.geometry as sg
import shapely.ops as so
import sys

def rotation_matrix(v1, v2):
    # formula for rotation matrix from:
    # https://math.stackexchange.com/questions/293116/rotating-one-3d-vector-to-another

    x = np.cross(v1,v2)
    x /= np.linalg.norm(x)

    a = np.array([
    [0, -1*x[2], x[1]],
    [x[2], 0, -1*x[0]],
    [-1*x[1], x[0], 0]])

    theta = np.dot(v1, v2)
    theta /= (np.linalg.norm(v1) * np.linalg.norm(v2) )
    theta = np.arccos(theta)

    return(linalg.expm(a * theta))

def align_n_to_c(atoms):
    com = np.mean(atoms, axis=0)
    #print("Center of mass",com)
    atoms_ = atoms - com
    v1 = atoms_[-1] - atoms_[0] # N to C terminus
    r1  = rotation_matrix(v1, np.array([0,1,0]))
    #print(r1)
    atoms_ = np.dot(atoms_, r1)
    r2 = rotation_matrix(np.array([atoms_[0,0],0,atoms_[0,2]]), np.array([0,0,1]))
    #print(r2)
    atoms_ = np.dot(atoms_, r2)
    return(atoms_ + com)

parser = PDBParser()
structure = parser.get_structure('PDB', sys.argv[1])
chain = structure[0]
coords = np.array([list(atom.get_vector()) for atom in chain.get_atoms()])

aligned_pts = align_n_to_c(coords)

out_prefix = sys.argv[2]

# create union of individual atoms
radius = 1.5 # approximate vdw radius in angstroms
space_filling = so.cascaded_union([sg.Point(i).buffer(radius) for i in aligned_pts])

# render and export
xs, ys = space_filling.simplify(0.3,preserve_topology=False).exterior.xy
fig, axs = plt.subplots()
axs.fill(xs, ys, alpha=0.5, fc='b', ec='k')
plt.axis('equal')
plt.savefig(out_prefix+'.pdf')
np.savetxt(out_prefix+'.csv',np.array([xs,ys]).T, delimiter=',')
