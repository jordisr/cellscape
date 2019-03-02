'''
pdb2svg:
    Scalable Vector Graphics for Macromolecular Structure

Author:
    Jordi Silvestre-Ryan (jordisr@berkeley.edu)

Usage of program and description of each command-line argument given with:
    python pdb2svg.py --help

Notes:
    In absence of view matrix from user, a view will be chosen by aligning the
    N-C terminal vector with the vertical axis.

'''

# Some color combinations I like:
# #FDB515' '#EE1F60' '#FDB515' '#3B7EA1'
# #F9E37E' '#E17272' '#F9E37E' '#F9977E'
# #276ab3' '#feb308' '#6fc276' '#ff9408'

from Bio.PDB import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import lines, text, cm
from matplotlib.colors import LinearSegmentedColormap
import shapely.geometry as sg
import shapely.ops as so
import os, sys, re, argparse, csv, pickle, colorsys, glob
from scipy import linalg

from parse_uniprot_xml import parse_xml

parser = argparse.ArgumentParser(description='Scalable Vector Graphics for Macromolecular Structure',  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# input pdb options
parser.add_argument('--pdb', help='Input PDB file')
parser.add_argument('--model', type=int, default=0, help='Model number in PDB to load')
parser.add_argument('--chain', default=['all'], help='Chain(s) in structure to outline', nargs='+')

# general input/output options
parser.add_argument('--view', help='File with output from PyMol get_view')
parser.add_argument('--uniprot', help='UniProt XML file to parse for sequence/domain/topology information')
parser.add_argument('--save', default='out', help='Prefix to save graphics')
parser.add_argument('--format', default='svg', help='Format to save graphics', choices=['svg','pdf','png'])
parser.add_argument('--export', action='store_true', help='Export Python object with structural information')
parser.add_argument('--look', help='Look in directory for structure .pdb, view matrix in .txt and UniProt .xml')

# visual style options
parser.add_argument('--residues', action='store_true', default=False, help='Draw residues separately and simulate surface rendering ')
parser.add_argument('--color_by', default='same',  choices=['same', 'chain', 'domain'], help='Color protein by chain')
parser.add_argument('--outline_domains', action='store_true', help='Outline each domain, implies --uniprot but not --residues')
parser.add_argument('--occlude', action='store_true', default=False, help='Occlude residues that are not visible')

# lower level graphics options
parser.add_argument('--radius', default=1.5, help='Space-filling radius, in angstroms', type=float)
parser.add_argument('--axes', action='store_true', default=False, help='Draw x and y axes around molecule')
parser.add_argument('--c', default=['#D3D3D3'], nargs='+', help='Set default color(s) in hex RGB')
parser.add_argument('--cmap', default='Set1', help='Set default color map')
parser.add_argument('--ec', default='k', help='Set default edge color')
parser.add_argument('--linewidth', default=0.7, type=float, help='Set default line width')

# residues to highlight separately
parser.add_argument('--highlight', type=int, help='List of residues to highlight',nargs='+')

# orientation and extra residues
parser.add_argument('--orientation', type=int, default=1, choices=[1,-1], help='Top-bottom orientation of protein (1:N>C or -1:C>N)')
parser.add_argument('--top-spacer', type=float, default=0, help='Placeholder at top of structure (length in nm)')
parser.add_argument('--bot-spacer', type=float, default=0, help='Placeholder at bottom of structure (length in nm)')

# arguments that need reworking
parser.add_argument('--topology', help='CSV-formatted file with topology boundaries')
parser.add_argument('--scale-bar', action='store_true', default=False, help='Draw a scale bar (experimental)')

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
    c3 = (h, min(l+(1-w)*l, 1), s)
    c2 = (h, l, s)
    c1 = (h, max(l-(1-w)*l, 0), s)
    # convert back to RGB and return colormap
    colors = [colorsys.hls_to_rgb(c[0], c[1], c[2]) for c in [c1, c2, c3]]
    return LinearSegmentedColormap.from_list(name, colors)

def rgb_to_cmap(h, w=0.3, name='test'):
    r = h[0]
    g = h[1]
    b = h[2]
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    # lighter and darker versions of color in HLS space
    c3 = (h, min(l+(1-w)*l, 1), s)
    c2 = (h, l, s)
    c1 = (h, max(l-(1-w)*l, 0), s)
    # convert back to RGB and return colormap
    colors = [colorsys.hls_to_rgb(c[0], c[1], c[2]) for c in [c1, c2, c3]]
    return LinearSegmentedColormap.from_list(name, colors)

def plot_polygon(poly, fc):
    axs = plt.gca()
    if isinstance(poly, sg.polygon.Polygon):
        xs, ys = poly.exterior.xy
        axs.fill(xs, ys, alpha=1, fc=fc, ec=args.ec, linewidth=args.linewidth, zorder=3)
    elif isinstance(poly, sg.multipolygon.MultiPolygon):
        for p in poly:
            xs, ys = p.exterior.xy
            axs.fill(xs, ys, alpha=1, fc=fc, ec=args.ec, linewidth=args.linewidth, zorder=3)
    return 0

def get_sequential_colors(n):
    if len(args.c) > 1:
        sequential_colors = args.c
    else:
        cmap = cm.get_cmap(args.cmap)
        sequential_colors = [cmap(x) for x in range(n)]
    return sequential_colors

def get_sequential_cmap(n):
    if len(args.c) > 1:
        cmap_list = [hex_to_cmap(c) for c in args.c]
    else:
        #cmap_x = np.linspace(0.0,1.0, len(chain_selection)) # continuous cmap
        #cmap_list = [rgb_to_cmap(cmap(x)) for x in cmap_x] # continuous cmap
        cmap = cm.get_cmap(args.cmap)
        cmap_list = [rgb_to_cmap(cmap(x)) for x in range(n)]
    return cmap_list

def get_residue_atoms(structure, residue):
    return np.array([list(a.get_vector()) for a in structure[r]])

def orientation_from_topology(topologies):
    first_ex_flag = True
    first_cy_flag = True
    first_he_flag = True

    for row in topologies:
        (description, start, end) = row

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
        orient_from_topo = 1
    elif first_ex[0] > first_cy[0]:
        orient_from_topo = -1

    return(orient_from_topo)

if __name__ == '__main__':

    if args.look:
        pdb_files = glob.glob(args.look+'/*.pdb')
        xml_files = glob.glob(args.look+'/*.xml')
        txt_files = glob.glob(args.look+'/*.txt')
        assert(len(pdb_files) > 0 or args.pdb)
        if len(xml_files) > 0:
            args.pdb = pdb_files[0]
        if len(xml_files) > 0:
            args.uniprot = xml_files[0]
        if len(txt_files) > 0:
            args.view = txt_files[0]

    # parse UniProt XML file if present
    if args.uniprot:
        up = parse_xml(args.uniprot)[0]
    elif args.color_by == 'domain' or args.outline_domains:
        sys.exit("Error: No domain information. Need to specify topology with UniProt XML file (--uniprot)")

    # infer orientation of protein from UniProt topology, if present
    if args.uniprot and len(up.topology) > 0:
        print(up.topology)
        orientation = orientation_from_topology(up.topology)
    else:
        orientation = args.orientation

    # open PDB structure
    parser = PDBParser()
    structure = parser.get_structure('PDB', args.pdb)
    model = structure[args.model]

    # select desired chains
    if args.chain[0] == 'all':
        chain_selection = [chain.id for chain in model.get_chains()]
    else:
        chain_selection = args.chain

    # apply rotation
    if args.view:
        # load view matrix
        view_mat = read_pymol_view(args.view)
        model.transform(view_mat,[0,0,0])
    else:
        # align N to C terminus
        untransformed_coords = np.concatenate([np.array([list(atom.get_vector()) for atom in model[chain].get_atoms()]) for chain in chain_selection])
        model.transform(align_n_to_c_mat(untransformed_coords, orientation),[0,0,0])

    # recenter coordinates on lower left edge of bounding box
    untransformed_coords = np.concatenate([np.array([list(atom.get_vector()) for atom in model[chain].get_atoms()]) for chain in chain_selection])
    offset_x = np.min(untransformed_coords[:,0])*1.01 # Shapely bug: non-noded intersection
    offset_y = np.min(untransformed_coords[:,1])
    model.transform(np.identity(3), -1*np.array([offset_x, offset_y, 0]))

    # calculate vertical offset for transmembrane proteins
    if args.uniprot and len(up.topology) > 0:
        topologies = up.topology
        tm_coordinates = []
        for t in topologies:
            (description, start, end) = t
            if description == 'Helical':
                for r in range(start, end+1):
                    tm_coordinates.append(get_residue_atoms(model[chain_selection[0]], r))
        tm_coordinates = np.concatenate(np.array(tm_coordinates))
        tm_com_y = np.mean(tm_coordinates[:,1])
        model.transform(np.identity(3), -1*np.array([0, tm_com_y+20, 0]))

    # global list of all atoms
    atom_coords = np.concatenate([np.array([list(atom.get_vector()) for atom in model[chain].get_atoms()]) for chain in chain_selection])

    # dict of dicts holding residue to np.array of atomic coordinates
    residue_to_atoms = dict()
    for chain in chain_selection:
        residue_to_atoms[chain] = dict()
        for residue in model[chain].get_residues():
            res_id = residue.get_full_id()[3][1]
            residue_to_atoms[chain][res_id] = np.array([list(r.get_vector()) for r in Selection.unfold_entities(residue,'A')])

    if args.uniprot and len(up.domains) > 0:

        # get mappings of residue to domain, empty if no UniProt data
        residue_to_domains = dict()
        d = 0
        prev_domain = 'None'
        for row in up.domain_segments:
            (name, start, end) = row
            if prev_domain != 'None':
                d += 1
            for r in range(start, end+1):
                residue_to_domains[r] = d
            prev_domain = name

        # retrieve domain residues and coordinates
        domain_atoms = []
        for row in up.domains:
            (name, start, end) = row
            this_domain = []
            for r in range(start, end+1):
                if r in residue_to_atoms[chain_selection[0]]:
                    this_domain.append(residue_to_atoms[chain_selection[0]][r])
                else:
                    print('WARNING: Missing residue',r,'in structure!')
            domain_atoms.append(np.concatenate(this_domain))

    # fire up a pyplot
    fig, axs = plt.subplots()
    axs.set_aspect('equal')

    if args.axes:
        axs.xaxis.grid(False)
        axs.yaxis.grid(True)
        axs.axes.xaxis.set_ticklabels([])
    else:
        plt.axis('off')
        plt.gca().set_axis_off()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
        plt.margins(0,0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())

    # unstructured regions, temporary fix of line at top or bottom of protein
    if args.orientation == 1:
        top_id = -1
        bot_id = 0
    else:
        top_id = 0
        bot_id = -1

    if args.top_spacer:
        top_spacer = lines.Line2D([atom_coords[top_id,0],atom_coords[top_id,0]], [atom_coords[top_id,1],atom_coords[top_id,1]+args.orientation*args.top_spacer*10], color=args.c[0], axes=axs, lw=10, zorder=0)
        axs.add_line(top_spacer)

    if args.bot_spacer:
        bot_spacer = lines.Line2D([atom_coords[bot_id,0],atom_coords[bot_id,0]], [atom_coords[bot_id,1],atom_coords[bot_id,1]-args.orientation*args.bot_spacer*10], color=args.c[0], axes=axs, lw=10, zorder=0)
        axs.add_line(bot_spacer)

    # main drawing routines
    if args.residues:
        atom_coords = np.array([])
        res_data = []
        for chain_i, chain in enumerate(chain_selection):
            for res_id, coords in residue_to_atoms[chain].items():
                if not args.uniprot:
                    domain_id = ""
                else:
                    domain_id = residue_to_domains.get(res_id, "")
                res_data.append((chain_i, res_id, coords, domain_id))
                if len(atom_coords) == 0:
                    atom_coords = coords
                else:
                    atom_coords = np.append(atom_coords,coords,axis=0)

        if args.color_by != 'same':
            if args.color_by == 'chain':
                color_on_col = 0
                cmap_list = get_sequential_cmap(len(chain_selection))
            elif args.color_by == 'domain':
                color_on_col = 3
                cmap_list = get_sequential_cmap(len(up.domains))
        else:
            cmap_list = [hex_to_cmap(args.c[0])]
            color_on_col = 0
        cmap_colors = len(cmap_list)

        z_coord = atom_coords[:,2]
        def rescale_coord(z):
            return (z-np.min(z_coord))/(np.max(z_coord)-np.min(z_coord))

        if args.occlude:
            # currently only works with outlining overlapping chains
            res_data_occlusion = []
            region_polygons = [[] for i in range(len(chain_selection))]
            view_object = None
            for row in reversed(sorted(res_data, key=lambda x: np.mean(x[2][:,2]))):
                coords = row[2]
                space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in coords])
                if not view_object:
                    view_object = space_filling
                else:
                    if view_object.disjoint(space_filling):
                        res_data_occlusion.append(row)
                        view_object = view_object.union(space_filling)
                        region_polygons[row[0]].append(space_filling)
                    elif view_object.contains(space_filling):
                        pass
                    else:
                        res_data_occlusion.append(row)
                        region_polygons[row[0]].append(space_filling.difference(view_object))
                        view_object = view_object.union(space_filling)

            sequential_colors = get_sequential_colors(len(chain_selection))
            for i, rp in enumerate(region_polygons):
                merged = so.cascaded_union(rp)
                plot_polygon(merged, fc=sequential_colors[i % len(sequential_colors)])

        else:
            for row in sorted(res_data, key=lambda x: np.mean(x[2][:,2]), reverse=False):
                coords = row[2]
                res_color = cmap_list[row[color_on_col] % cmap_colors](rescale_coord(np.mean(coords[:,2])))
                space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in coords])
                plot_polygon(space_filling, fc=res_color)

    else:
        if args.outline_domains:
            sequential_colors = get_sequential_colors(len(up.domains))
            for i,coords in enumerate(domain_atoms):
                domain_coords = coords
                space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in domain_coords])
                plot_polygon(space_filling, fc=sequential_colors[i % len(sequential_colors)])

        else:
            space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in atom_coords])
            plot_polygon(space_filling, fc=args.c[0])

    if args.highlight:
        # draws those residues separately on top of previous polygons
        highlight_res = [residue_to_atoms[int(i)] for i in args.highlight]
        for v in highlight_res:
            res_coords = v
            space_filling = so.cascaded_union([sg.Point(i).buffer(args.radius) for i in res_coords])
            plot_polygon(space_filling, fc='r')

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
    plt.savefig(out_prefix+'.'+args.format, transparent=True, pad_inches=0, bbox_inches='tight', dpi=300)

    # output summary line
    image_width = np.max(atom_coords[:,0]) - np.min(atom_coords[:,0])
    image_height = np.max(atom_coords[:,1]) - np.min(atom_coords[:,1])
    print('\t'.join(map(str, [args.pdb, args.save+'.'+args.format, image_height, image_width])))

    if args.export:
        data = {'polygons':[], 'width':image_width, 'height':image_height}
        import matplotlib
        for o in plt.gca().findobj(matplotlib.patches.Polygon):
            data['polygons'].append(o)
            print(o)
        with open(args.save+'.pickle','wb') as f:
            pickle.dump(data, f)
