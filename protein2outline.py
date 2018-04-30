'''
Generate 2D space-filling outline from a PDB structure.

In absence of view matrix, a view will be chosen by aligning the N-C terminal
vector with the vertical axis.

'''

from Bio.PDB import *
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import shapely.geometry as sg
import shapely.ops as so
import os, sys, re, argparse

# parse command line arguments
parser = argparse.ArgumentParser(description='Produce space-filling vector graphics outline of PDB structure')
parser.add_argument('--pdb', help='PDB structure to load')
parser.add_argument('--view', help='File with view matrix')
parser.add_argument('--model', default=0, help='Model in PDB to load')
parser.add_argument('--chain', default='A', help='Chain in PDB to outline')
parser.add_argument('--recenter', action='store_true', default=False, help='Recenter atomic coordinates')
parser.add_argument('--highlight', type=int, help='Residues to highlight',nargs='+')
parser.add_argument('--domains', type=int, help='List of domain boundaries', nargs='+')
#args = parser.parse_args()

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

def rotate_atoms(atoms):
    test_mat = np.array([
    [-0.350587279, 0.472307622, -0.808711231],
    [0.560661256, 0.797530472, 0.222723737],
    [0.750165045, -0.375327766, -0.544408023]
    ])
    atoms_ = np.dot(atoms, test_mat)
    return(atoms_)

if __name__ == '__main__':

    parser = PDBParser()
    structure = parser.get_structure('PDB', sys.argv[1])
    chain = structure[0]
    coords = np.array([list(atom.get_vector()) for atom in chain.get_atoms()])

    for residue in chain.get_residues():
        res_atoms = residue.get_unpacked_list()
        for atom in res_atoms:
            print(atom.get_vector())

    #aligned_pts = align_n_to_c(coords)
    aligned_pts = rotate_atoms(coords)

    out_prefix = sys.argv[2]

    # create union of individual atoms
    radius = 1.5 # approximate vdw radius in angstroms
    space_filling = so.cascaded_union([sg.Point(i).buffer(radius) for i in aligned_pts])

    # render and export
    xs, ys = space_filling.simplify(0.3,preserve_topology=False).exterior.xy
    fig, axs = plt.subplots()
    axs.fill(xs, ys, alpha=0.5, fc='b', ec='k')
    plt.axis('equal')
    plt.axis('off')
    plt.savefig(out_prefix+'.pdf',transparent=True)
    np.savetxt(out_prefix+'.csv',np.array([xs,ys]).T, delimiter=',')
