'''
pdb2svg:
    Scalable Vector Graphics for Macromolecular Structure

Author:
    Jordi Silvestre-Ryan (jordisr@berkeley.edu)

Usage of program and description of each command-line argument given with:
    python pdb2svg.py --help

Notes:
    In absence of view matrix from user, a view will be chosen by aligning the
    N-C terminal vector with the vertical axis (in progress).

TO DO:
- work on backbone ribbon option with splines
- represent unstructured regions
- write tutorial/documentation
- decide where to put style options
- compatibility with other molecular graphics programs (Chimera, VMD?)
'''

from Bio.PDB import *
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from matplotlib import lines, text
import shapely.geometry as sg
import shapely.ops as so
import os, sys, re, argparse, csv

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

# visual options
parser.add_argument('--scale-bar', action='store_true', default=False, help='Draw a scale bar')
parser.add_argument('--axes', action='store_true', default=False, help='Draw x and y axes around molecule')

# draw separate polygon around each residue, entire protein, or each domain
parser.add_argument('--all', action='store_true', default=False, help='Draw all residues separately (overrides --domains and --backbone)')
parser.add_argument('--domains', help='CSV-formatted file with region/domain boundaries')
parser.add_argument('--outline', action='store_true', default=True, help='Draw one outline for entire structure (default behavior)')

# experimental arguments, override other options
parser.add_argument('--backbone', default=False, choices=['all','ca'], help='(Experimental) Draw backbone with splines')
parser.add_argument('--unstructured', action='store_true', default=False, help='(Experimental) Extra regions')

args = parser.parse_args()

# arguments to be incorporated later
parser.add_argument('--model', default=0, help='Model number in PDB to load')
parser.add_argument('--chain', default='A', help='Chain(s) in structure to outline', nargs='+')
parser.add_argument('--recenter', action='store_true', default=False, help='Recenter atomic coordinates')

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

    if args.backbone:
        if args.backbone == 'all':
            atom_group = ('N','CA','C', 'O')
        elif args.backbone == 'ca':
            atom_group = ('CA')
        backbone_atoms = []
        for atom in chain.get_atoms():
            atom_id = atom.get_full_id()[-1][0]
            if atom_id in atom_group:
                backbone_atoms.append(list(atom.get_vector()))
        backbone_atoms = np.array(backbone_atoms)

    # dictionary holding residue to np.array of atomic coordinates
    residue_to_atoms = dict()
    for residue in chain.get_residues():
        res_id = residue.get_full_id()[3][1]
        residue_to_atoms[res_id] = np.array([list(r.get_vector()) for r in Selection.unfold_entities(residue,'A')])

    if args.domains:
        domain_atoms = []
        with open(args.domains) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                (start,end) = (int(row['res_start']),int(row['res_end']))
                domain_atoms.append(np.concatenate([residue_to_atoms[r] for r in range(start,end)]))

    mat = read_pymol_view(args.view)
    aligned_pts = np.dot(coords,mat)
    #aligned_pts = align_n_to_c(coords)

    # fire up a pyplot
    fig, axs = plt.subplots()
    #plt.axis('equal')
    axs.set_aspect('equal')

    if args.axes:
        axs.xaxis.grid(False)
        axs.yaxis.grid(True)
        axs.axes.xaxis.set_ticklabels([])
    else:
        plt.axis('off')

    # create space filling representation
    if args.domains:
        sequential_colors = ['#FDB515','#00356B']
        for i,coords in enumerate(domain_atoms):
            domain_coords = np.dot(coords,mat)
            space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in domain_coords])
            xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
            axs.fill(xs, ys, alpha=1, fc=sequential_colors[i % len(sequential_colors)], ec='k')
    elif args.all:
        for i,coords in residue_to_atoms.items():
            domain_coords = np.dot(coords,mat)
            space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in domain_coords])
            xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
            axs.fill(xs, ys, alpha=1, fc='#D3D3D3', ec='#A9A9A9')
    elif args.outline:
        space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in aligned_pts])
        xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
        axs.fill(xs, ys, alpha=1, fc='#377eb8', ec='k')

    if args.backbone:
        # backbone rendering, needs work
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

    if args.highlight:
        # draws those residues separately on top of previous polygons
        #highlight_com = np.dot(np.array(list(highlight_res.values())),mat)
        #plt.scatter(highlight_com[:,0],highlight_com[:,1], c='k')
        highlight_res = [residue_to_atoms[int(i)] for i in args.highlight]
        for v in highlight_res:
            res_coords = np.dot(v,mat)
            space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in res_coords])
            xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
            axs.fill(xs, ys, alpha=1, fc='r', ec='k')

    if args.scale_bar:
        bar_length = 1*10
        bar_pos_x = aligned_pts[0,0]
        bar_pos_y = aligned_pts[0,1]
        scale_bar = lines.Line2D([bar_pos_x,bar_pos_x], [bar_pos_y,bar_pos_y+bar_length], color='black', axes=axs, lw=5)
        axs.add_line(scale_bar)
        # legend for scale bar
        #legend = text.Text(0,100, 'text label', ha='left', va='bottom', axes=axs)
        #axs.add_artist(legend)

    # output coordinates and vector graphics
    out_prefix = args.save
    plt.savefig(out_prefix+'.'+args.format,transparent=True)
    #np.savetxt(out_prefix+'.csv',np.array([xs,ys]).T, delimiter=',') # save coordinates in csv
