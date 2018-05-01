'''
pdb2svg: Scalable Vector Graphics for Macromolecular Structure

molecule2outline, pdb2outline, outline-molecule, MolecularOutline

Author: Jordi Silvestre-Ryan (jordisr@berkeley.edu)

N.B. In absence of view matrix from user, a view will be chosen by aligning the N-C terminal
vector with the vertical axis.

TO DO:
- reorganize code
- read domains/colors from file?
- add backbone ribbon option (with splines?)
- compatibility with other molecular graphics programs (Chimera, VMD?)
- write tutorial/documentation
'''

from Bio.PDB import *
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import shapely.geometry as sg
import shapely.ops as so
import os, sys, re, argparse

def check_simplify(value):
    fvalue = float(value)
    if fvalue > 0 and fvalue < 2:
        return True
    else:
        raise argparse.ArgumentTypeError("Simplify called with value of %s. Choose a value between 0-2." % value)

# parse command line arguments
parser = argparse.ArgumentParser(description='Scalable Vector Graphics for Macromolecular Structure',  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--pdb', help='Input PDB file', required=True)
parser.add_argument('--view', help='File with PyMol view', required=True)
parser.add_argument('--save', default='out', help='Prefix to save graphics')
parser.add_argument('--radius', default=1.5, help='Space-filling radius, in angstroms', type=float)
parser.add_argument('--simplify', default=0, help='Amount to simplify resulting polygons', type=check_simplify)
parser.add_argument('--highlight', type=int, help='List of residues to highlight',nargs='+')
parser.add_argument('--format', default='svg', help='Format to save graphics', choices=['svg','pdf'])

# experimental arguments
parser.add_argument('--test', action='store_true', default=False, help='Experimental!')
parser.add_argument('--all_atom', action='store_true', default=False, help='Experimental! Draw each residue')
parser.add_argument('--backbone', action='store_true', default=False, help='Experimental!')
parser.add_argument('--ca_backbone', action='store_true', default=False, help='Experimental!')


args = parser.parse_args()

# arguments to be incorporated later
parser.add_argument('--domains', type=int, help='List of domain boundaries, e.g. 1 10 11 20', nargs='+')
parser.add_argument('--model', default=0, help='Model number in PDB to load')
parser.add_argument('--chain', default='A', help='Chain(s) in structure to outline', nargs='+')
parser.add_argument('--recenter', action='store_true', default=False, help='Recenter atomic coordinates')

def rotation_matrix(v1, v2):
    # formula for rotation matrix from:McTest
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

# read rotation matrix from PyMol get_view command
def read_pymol_view(file):
    matrix = []
    with open(file,'r') as view:
        for line in view:
            fields =line.split(',')
            if len(fields) == 4:
                matrix.append(list(map(float,fields[:3])))
    return(np.array(matrix)[:3])

if __name__ == '__main__':

    parser = PDBParser()
    structure = parser.get_structure('PDB', args.pdb)
    chain = structure[0]
    coords = np.array([list(atom.get_vector()) for atom in chain.get_atoms()])

    # set up residue highlights
    highlight_res = dict()

    # backbone atoms
    backbone_atoms = []
    ca_atoms = []
    for atom in chain.get_atoms():
        atom_id = atom.get_full_id()[-1][0]
        if atom_id in ('N','CA','C', 'O'):
            backbone_atoms.append(list(atom.get_vector()))
        if atom_id == 'CA':
            ca_atoms.append(list(atom.get_vector()))
    backbone_atoms = np.array(backbone_atoms)
    ca_atoms = np.array(ca_atoms)

    # collect highlighted residues
    if args.highlight:
        for residue in chain.get_residues():
            res_id = residue.get_full_id()[3][1]
            if res_id in args.highlight:
                #highlight_res[res_id] = np.mean(np.array([list(r.get_vector()) for r in Selection.unfold_entities(residue,'A')]), axis=0)
                highlight_res[res_id] = np.array([list(r.get_vector()) for r in Selection.unfold_entities(residue,'A')])
                #highlight_res[res_id] = np.array(list(residue['CA'].get_vector()))

    if args.test:
        test_domains = [(35,144), (146,237), (238, 322), (324, 415), (416, 498), (502, 593), (594,677)] # ceacam5 domain boundaries
        residue_to_atoms = dict()
        for residue in chain.get_residues():
            res_id = residue.get_full_id()[3][1]
            residue_to_atoms[res_id] = np.array([list(r.get_vector()) for r in Selection.unfold_entities(residue,'A')])
        domain_atoms = []
        for start,end in test_domains:
            domain_atoms.append(np.concatenate([residue_to_atoms[r] for r in range(start,end)]))

    # transform and project atoms
    #test_mat = np.array([
    #    [-0.350587279, 0.472307622, -0.808711231],
    #    [0.560661256, 0.797530472, 0.222723737],
    #    [0.750165045, -0.375327766, -0.544408023]
    #])
    mat = read_pymol_view(args.view)
    aligned_pts = np.dot(coords,mat)
    #aligned_pts = align_n_to_c(coords)

    # fire up a pyplot
    fig, axs = plt.subplots()
    plt.axis('equal')
    plt.axis('off')

    # create space filling representation
    if not args.test and not args.all_atom:
        space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in aligned_pts])
        xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
        axs.fill(xs, ys, alpha=1, fc='#377eb8', ec='k')
    elif args.test and not args.all_atom:
        #sequential_colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33']
        sequential_colors = ['#FDB515','#00356B']
        for i,coords in enumerate(domain_atoms):
            domain_coords = np.dot(coords,mat)
            space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in domain_coords])
            xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
            axs.fill(xs, ys, alpha=1, fc=sequential_colors[i % len(sequential_colors)], ec='k')
    elif args.all_atom:
        for i,coords in residue_to_atoms.items():
            domain_coords = np.dot(coords,mat)
            space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in domain_coords])
            xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
            axs.fill(xs, ys, alpha=1, fc='#D3D3D3', ec='#A9A9A9')

    if args.backbone and not args.ca_backbone:
        from scipy import interpolate
        from scipy.signal import savgol_filter
        backbone_atoms = np.dot(backbone_atoms, mat)
        fx = savgol_filter(backbone_atoms[:,0],21,3)
        fy = savgol_filter(backbone_atoms[:,1],21,3)
        plt.plot(fx, fy, c='k')
        # now do optional interpolation for final coordinates
        #tck, u = interpolate.splprep([fx, fy],s=3)
        #unew = np.arange(0, 1.01, 0.01)
        #out = interpolate.splev(unew, tck)
        #plt.plot(out[0], out[1])
    elif args.ca_backbone and not args.backbone:
        from scipy import interpolate
        from scipy.signal import savgol_filter
        ca_atoms = np.dot(ca_atoms, mat)
        plt.plot(ca_atoms[:,0], ca_atoms[:,1], c='k')
        #fx = savgol_filter(ca_atoms[:,0],5,3)
        #fy = savgol_filter(ca_atoms[:,1],5,3)
        #plt.plot(fx, fy, c='k')
        # now do optional interpolation for final coordinates
        #tck, u = interpolate.splprep([fx,fy],s=0)
        #unew = np.arange(0, 1.01, 0.01)
        #out = interpolate.splev(unew, tck)
        #plt.plot(out[0], out[1], c='k')

    if args.highlight:
        #highlight_com = np.dot(np.array(list(highlight_res.values())),mat)
        #plt.scatter(highlight_com[:,0],highlight_com[:,1], c='k')
        for k,v in highlight_res.items():
            res_coords = np.dot(v,mat)
            space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in res_coords])
            xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
            axs.fill(xs, ys, alpha=1, fc='r', ec='k')

    # output coordinates and vector graphics
    out_prefix = args.save
    plt.savefig(out_prefix+'.'+args.format,transparent=True)
    #np.savetxt(out_prefix+'.csv',np.array([xs,ys]).T, delimiter=',') # save coordinates in csv
