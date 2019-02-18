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
from matplotlib.colors import LinearSegmentedColormap
import colorsys
import shapely.geometry as sg
import shapely.ops as so
import os, sys, re, argparse, csv, pickle
from scipy.signal import savgol_filter
from scipy import interpolate

def check_simplify(value):
    fvalue = float(value)
    if fvalue > 0 and fvalue < 2:
        return True
    else:
        raise argparse.ArgumentTypeError("Simplify called with value of %s. Choose a value between 0-2." % value)

parser = argparse.ArgumentParser(description='Scalable Vector Graphics for Macromolecular Structure',  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# input pdb options
parser.add_argument('--pdb', help='Input PDB file', required=True)
parser.add_argument('--model', type=int, default=0, help='Model number in PDB to load')
parser.add_argument('--chain', default=['A'], help='Chain(s) in structure to outline', nargs='+')

# general input/output options
parser.add_argument('--view', help='File with output from PyMol get_view')
parser.add_argument('--save', default='out', help='Prefix to save graphics')
parser.add_argument('--format', default='svg', help='Format to save graphics', choices=['svg','pdf'])
parser.add_argument('--export', action='store_true', help='Export Python object with structural information')

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
parser.add_argument('--depth', action='store_true', default=False, help='Experimental rendering')
parser.add_argument('--outline', action='store_true', default=True, help='Draw one outline for entire structure (default behavior)')

# orientation and extra residues
parser.add_argument('--orientation', type=int, default=1, choices=[1,-1], help='Top-bottom orientation of protein (1:N>C or -1:C>N)')
parser.add_argument('--top-spacer', type=float, default=0, help='Placeholder at top of structure (length in nm)')
parser.add_argument('--bot-spacer', type=float, default=0, help='Placeholder at bottom of structure (length in nm)')

# arguments that need reworking
parser.add_argument('--domains', help='CSV-formatted file with region/domain boundaries')
parser.add_argument('--recenter', type=int, default=0, help='Recenter atomic coordinates on this residue')
parser.add_argument('--topology', help='CSV-formatted file with topology boundaries')

args = parser.parse_args()

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

def align_n_to_c_mat(atoms, sign=1):
    com = np.mean(atoms, axis=0)
    atoms_ = atoms - com
    v1 = sign*atoms_[-1] - sign*atoms_[0] # N to C terminus
    r1  = rotation_matrix(v1, np.array([0,1,0]))
    return(r1)

# read rotation matrix from PyMol get_view command
def read_pymol_view(file):
    matrix = []
    with open(file,'r') as view:
        for line in view:
            fields =line.split(',')
            if len(fields) == 4:
                matrix.append(list(map(float,fields[:3])))
    return(np.array(matrix)[:3])

# define custom Lighter Color => Color => Darker Color cmap
def hex_to_cmap(h, w=0.3, name='test'):
    r = int(h[1:3], 16)
    g = int(h[3:5], 16)
    b = int(h[5:7], 16)
    h, l, s = colorsys.rgb_to_hls(r/255,g/255,b/255)
    # lighter and darker versions of color in HLS space
    c1 = (h, min(l+(1-w)*l, 1), s)
    c2 = (h, l, s)
    c3 = (h, max(l-(1-w)*l, 0), s)
    # convert back to RGB and return colormap
    colors = [colorsys.hls_to_rgb(c[0], c[1], c[2]) for c in [c1, c2, c3]]
    return LinearSegmentedColormap.from_list(name, colors)

if __name__ == '__main__':

    parser = PDBParser()
    structure = parser.get_structure('PDB', args.pdb)
    model = structure[args.model]

    # select desired chains
    if args.chain[0] == 'all':
        chain_selection = [chain.id for chain in model.get_chains()]
    else:
        chain_selection = args.chain

    # set up residue highlights
    highlight_res = dict()

    if args.view:
        # load view matrix
        view_mat = read_pymol_view(args.view)
        model.transform(view_mat,[0,0,0])
    else:
        # align N to C terminus
        model.transform(align_n_to_c_mat(atom_coords,orient_from_topo),[0,0,0])

    # rewrite so this line isn't in twice
    atom_coords = np.concatenate([np.array([list(atom.get_vector()) for atom in model[chain].get_atoms()]) for chain in chain_selection])

    # recenter coordinates on residue (useful for orienting transmembrane proteins)
    if args.recenter:
        offset_res_id = args.recenter
        # accessing chain by residue id seems to be problematic?
        for residue in model.get_residues():
            res_id = residue.get_full_id()[3][1]
            if res_id == offset_res_id:
                (offset_x, offset_y, _) = np.mean(np.array([list(r.get_vector()) for r in Selection.unfold_entities(residue,'A')]),axis=0)
    else:
        offset_x = atom_coords[0,0]
        offset_y = atom_coords[0,1]
        model.transform(np.identity(3), [-1*offset_x,-1*offset_y,0])

    #atom_coords = np.array([list(atom.get_vector()) for atom in model.get_atoms()])
    atom_coords = np.concatenate([np.array([list(atom.get_vector()) for atom in model[chain].get_atoms()]) for chain in chain_selection])

    # dict of dicts holding residue to np.array of atomic coordinates
    residue_to_atoms = dict()
    for chain in chain_selection:
        residue_to_atoms[chain] = dict()
        for residue in model[chain].get_residues():
            res_id = residue.get_full_id()[3][1]
            #print(residue.get_full_id())
            residue_to_atoms[chain][res_id] = np.array([list(r.get_vector()) for r in Selection.unfold_entities(residue,'A')])

    if args.domains:
        domain_atoms = []
        with open(args.domains) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                (start,end) = (int(row['res_start']),int(row['res_end']))
                this_domain = []
                for r in range(start,end):
                    if r in residue_to_atoms:
                        this_domain.append(residue_to_atoms[r])
                    else:
                        print('WARNING: Missing residue',r,'in structure!')
                domain_atoms.append(np.concatenate(this_domain))

    if args.topology:
        # not currently used for automatic orienting
        topologies = []
        first_ex_flag = True
        first_cy_flag = True
        first_he_flag = True
        with open(args.topology) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                (start,end,description) = (int(row['res_start']),int(row['res_end']),row['description'])
                topologies.append((start,end,description))
                if description == 'Extracellular' and first_ex_flag:
                    first_ex = (start, end)
                    first_ex_flag = False
                elif description == 'Helical' and first_he_flag:
                    first_he = (start, end)
                    first_he_flag = False
                elif description == 'Cytoplasmic' and first_cy_flag:
                    first_cy = (start, end)
                    first_cy_flag = False
        if first_ex[0] < first_cy[0]:
            orient_from_topo = -1
        elif first_ex[0] > first_cy[0]:
            orient_from_topo = 1
        if orient_from_topo == 1:
            tm_start = first_cy[1]
        elif orient_from_topo == -1:
            tm_start = first_cy[0]

        def safe_bounds(n, residue_to_atoms):
            res_in_struct = list(residue_to_atoms.keys())
            if n < np.min(res_in_struct):
                return(np.min(res_in_struct))
            elif n > np.max(res_in_struct):
                return(np.max(res_in_struct))

        offset_x = residue_to_atoms[safe_bounds(tm_start, residue_to_atoms)][0][0]
        offset_y = residue_to_atoms[safe_bounds(tm_start, residue_to_atoms)][0][1]

        chain.transform(np.identity(3), [-1*offset_x,-1*offset_y,0])
        atom_coords = np.array([list(atom.get_vector()) for atom in chain.get_atoms()])

    # fire up a pyplot
    fig, axs = plt.subplots()
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
        #sequential_colors = ['#FDB515','#00356B'] # berkeley
        #sequential_colors = ['#B22234','#FFFFFF','#3C3B6E']
        for i,coords in enumerate(domain_atoms):
            #domain_coords = np.dot(coords,mat)
            domain_coords = coords
            space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in domain_coords])
            xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
            #fx = savgol_filter(xs,11,2)
            #fy = savgol_filter(ys,11,2)
            #tck, u = interpolate.splprep([fx, fy],s=3)
            #unew = np.arange(0, 1.01, 0.01)
            #out = interpolate.splev(unew, tck)
            axs.fill(xs, ys, alpha=1, fc=sequential_colors[i % len(sequential_colors)], ec='k',zorder=2)

    elif args.depth:

        # only use atoms that will be drawn for outline (which ones are not here?)
        atom_coords = np.array([])
        res_data = []
        for chain_i, chain in enumerate(chain_selection):
            for res_id, coords in residue_to_atoms[chain].items():
                res_data.append((chain_i, res_id, coords))
                if len(atom_coords) == 0:
                    atom_coords = coords
                else:
                    atom_coords = np.append(atom_coords,coords,axis=0)

        # draw outline in the back first
        space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in atom_coords])
        try:
            xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
            axs.fill(xs, ys, alpha=1, fc='w', ec='k', zorder=1)
        except:
            pass

        # color maps for different chains
        #cmap_names = ['Blues', 'Oranges', 'Greens', 'Reds', 'Purples']
        #cmap_list = [cm.get_cmap(x) for x in cmap_names]
        cmap_list = [hex_to_cmap(c) for c in ['#276ab3', '#feb308', '#6fc276', '#ff9408']]
        cmap_colors = len(cmap_list)

        z_coord = atom_coords[:,2]
        def rescale_coord(z):
            return (z-np.min(z_coord))/(np.max(z_coord)-np.min(z_coord))

        for row in sorted(res_data, key=lambda x: np.mean(x[2][:,2]), reverse=True):
            coords = row[2]
            res_color = cmap_list[row[0] % cmap_colors](rescale_coord(np.mean(coords[:,2])))
            space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in coords])
            xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
            axs.fill(xs, ys, alpha=1, fc=res_color, ec='#202020', zorder=2,linewidth=0.4)

        #print(np.min(atom_coords[:,0]), np.max(atom_coords[:,0]), np.min(atom_coords[:,1]), np.max(atom_coords[:,1]))

    elif args.outline:
        space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in atom_coords])
        xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
        axs.fill(xs, ys, alpha=1, fc=args.c, ec='k', zorder=1)

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
        bar_length = 10
        bar_pos_x = np.min(atom_coords,axis=0)[0]
        bar_pos_y = atom_coords[0,1]
        scale_bar = lines.Line2D([bar_pos_x,bar_pos_x], [bar_pos_y,bar_pos_y+bar_length], color='black', axes=axs, lw=5)
        axs.add_line(scale_bar)
        # legend for scale bar
        #legend = text.Text(bar_pos_x,bar_pos_y+bar_length, '1 nm', ha='left', va='bottom', axes=axs)
        #axs.add_artist(legend)

    # output coordinates and vector graphics
    out_prefix = args.save
    plt.gca().set_axis_off()
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
    plt.margins(0,0)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.savefig(out_prefix+'.'+args.format, transparent=True, pad_inches=0, bbox_inches='tight')

    if args.export:
        # save 2d paths and protein metadata to a python pickle object0
        space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in atom_coords])
        xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
        outline = np.array([xs,ys])

        height = np.max(atom_coords[:,1]) - np.min(atom_coords[:,1])
        width = np.max(atom_coords[:,0]) - np.min(atom_coords[:,0])

        data = {
        'name': args.save,
        'outline': outline,
        'center': args.recenter,
        'height': height,
        'width':width,
        'orientation':args.orientation
        }

        if args.domains:
            domain_paths = []
            for i,coords in enumerate(domain_atoms):
                domain_coords = coords
                space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in domain_coords])
                xs, ys = space_filling.simplify(args.simplify,preserve_topology=False).exterior.xy
                domain_paths.append((xs,ys))
            data['domain_paths'] = domain_paths

        with open(args.save+'.pickle','wb') as f:
            pickle.dump(data, f)
