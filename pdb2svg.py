'''
pdb2svg:
    Scalable Vector Graphics for Macromolecular Structure

Author:
    Jordi Silvestre-Ryan (jordisr@berkeley.edu)

Usage of program and description of each command-line argument given with:
    python pdb2svg.py --help

Notes:
    In absence of view matrix from user, a view will be chosen by aligning the
    N-C terminal vector with the vertical axis (feature in proress).

'''

from Bio.PDB import *
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from matplotlib import lines, text, cm
import shapely.geometry as sg
import shapely.ops as so
import os, sys, re, argparse, csv

def check_simplify(value):
    fvalue = float(value)
    if fvalue > 0 and fvalue < 2:
        return True
    else:
        raise argparse.ArgumentTypeError("Simplify called with value of %s. Choose a value between 0-2." % value)

parser = argparse.ArgumentParser(description='Scalable Vector Graphics for Macromolecular Structure',  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# input/output options
parser.add_argument('--pdb', help='Input PDB file', required=True)
parser.add_argument('--view', help='File with PyMol view', required=True)
parser.add_argument('--save', default='out', help='Prefix to save graphics')
parser.add_argument('--format', default='svg', help='Format to save graphics', choices=['svg','pdf'])

# visual style options
parser.add_argument('--radius', default=1.5, help='Space-filling radius, in angstroms', type=float)
parser.add_argument('--simplify', default=0, help='Amount to simplify resulting polygons', type=check_simplify)
parser.add_argument('--scale-bar', action='store_true', default=False, help='Draw a scale bar')
parser.add_argument('--axes', action='store_true', default=False, help='Draw x and y axes around molecule')
parser.add_argument('--c', default='#D3D3D3', help='Color (if used)')
parser.add_argument('--cmap', default='jet', help='Colormap (if used)')

# residues to highlight separately
parser.add_argument('--highlight', type=int, help='List of residues to highlight',nargs='+')

# draw separate polygon around each residue, entire protein, or each domain
parser.add_argument('--all', action='store_true', default=False, help='Draw all residues separately (overrides --domains and --backbone)')
parser.add_argument('--domains', help='CSV-formatted file with region/domain boundaries')
parser.add_argument('--outline', action='store_true', default=True, help='Draw one outline for entire structure (default behavior)')

# experimental arguments, override other options
parser.add_argument('--backbone', default=False, choices=['all','ca'], help='(Experimental) Draw backbone with splines')

# orientation and extra residues
parser.add_argument('--orientation', type=int, default=1, choices=[1,-1], help='Top-bottom orientation of protein (1:N>C or -1:C>N)')
parser.add_argument('--top-spacer', type=float, default=0, help='Placeholder at top of structure (length in nm)')
parser.add_argument('--bot-spacer', type=float, default=0, help='Placeholder at bottom of structure (length in nm)')

parser.add_argument('--recenter', type=int, default=0, help='Recenter atomic coordinates on this residue')

args = parser.parse_args()

# arguments to be incorporated later
parser.add_argument('--model', default=0, help='Model number in PDB to load')
parser.add_argument('--chain', default='A', help='Chain(s) in structure to outline', nargs='+')

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

    # set up residue highlights
    highlight_res = dict()

    # load view matrix
    view_mat = read_pymol_view(args.view)

    # rotate coordinates with view matrix
    chain.transform(view_mat,[0,0,0])
    atom_coords = np.array([list(atom.get_vector()) for atom in chain.get_atoms()])

    # recenter coordinates on residue (useful for orienting transmembrane proteins)
    if args.recenter:
        offset_res_id = args.recenter
        # accessing chain by residue id seems to be problematic?
        for residue in chain.get_residues():
            res_id = residue.get_full_id()[3][1]
            if res_id == offset_res_id:
                (offset_x, offset_y, _) = np.mean(np.array([list(r.get_vector()) for r in Selection.unfold_entities(residue,'A')]),axis=0)
    else:
        offset_x = atom_coords[0,0]
        offset_y = atom_coords[0,1]
    chain.transform(np.identity(3), [-1*offset_x,-1*offset_y,0])
    atom_coords = np.array([list(atom.get_vector()) for atom in chain.get_atoms()])

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

    # unstructured regions, temporary fix of line at top or bottom of protein
    if args.orientation == 1:
        top_id = -1
        bot_id = 0
    else:
        top_id = 0
        bot_id = -1

    if args.top_spacer:
        top_spacer = lines.Line2D([atom_coords[top_id,0],atom_coords[top_id,0]], [atom_coords[top_id,1],atom_coords[top_id,1]+args.orientation*args.top_spacer*10], color=args.c, axes=axs, lw=10, zorder=0)
        axs.add_line(top_spacer)

    if args.bot_spacer:
        bot_spacer = lines.Line2D([atom_coords[bot_id,0],atom_coords[bot_id,0]], [atom_coords[bot_id,1],atom_coords[bot_id,1]-args.orientation*args.bot_spacer*10], color=args.c, axes=axs, lw=10, zorder=0)
        axs.add_line(bot_spacer)

    # create space filling representation
    if args.domains:
        # set color scheme
        cmap = cm.get_cmap(args.cmap)
        cmap_x = np.linspace(0.0,1.0,len(domain_atoms))
        sequential_colors = [cmap(x) for x in cmap_x]
        #sequential_colors = ['#FDB515','#00356B']
        for i,coords in enumerate(domain_atoms):
            #domain_coords = np.dot(coords,mat)
            domain_coords = coords
            space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in domain_coords])
            xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
            axs.fill(xs, ys, alpha=1, fc=sequential_colors[i % len(sequential_colors)], ec='k')
    elif args.all:
        for i,coords in residue_to_atoms.items():
            #domain_coords = np.dot(coords,mat)
            domain_coords = coords
            space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in domain_coords])
            xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
            axs.fill(xs, ys, alpha=1, fc='#D3D3D3', ec='#A9A9A9')
    elif args.outline:
        space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in atom_coords])
        xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
        axs.fill(xs, ys, alpha=1, fc=args.c, ec='k')

    if args.backbone:
        # backbone rendering, needs work
        from scipy import interpolate
        from scipy.signal import savgol_filter
        #backbone_atoms = np.dot(backbone_atoms, mat)
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
            #res_coords = np.dot(v,mat)
            res_coords = v
            space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in res_coords])
            xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
            axs.fill(xs, ys, alpha=1, fc='r', ec='k')

    if args.scale_bar:
        bar_length = 1*10
        bar_pos_x = np.min(atom_coords,axis=0)[0]
        bar_pos_y = atom_coords[0,1]
        scale_bar = lines.Line2D([bar_pos_x,bar_pos_x], [bar_pos_y,bar_pos_y+bar_length], color='black', axes=axs, lw=5)
        axs.add_line(scale_bar)
        # legend for scale bar
        #legend = text.Text(bar_pos_x,bar_pos_y+bar_length, '1 nm', ha='left', va='bottom', axes=axs)
        #axs.add_artist(legend)

    # output coordinates and vector graphics
    out_prefix = args.save
    plt.savefig(out_prefix+'.'+args.format,transparent=True)
    #np.savetxt(out_prefix+'.csv',np.array([xs,ys]).T, delimiter=',') # save coordinates in csv
