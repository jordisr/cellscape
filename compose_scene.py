import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import lines, text, cm
import shapely.geometry as sg
import shapely.ops as so
import os, sys, argparse, pickle
import glob
import csv

parser = argparse.ArgumentParser(description='Structure SVG compositing',  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--files', nargs='+', help='Pickled objects to load')
parser.add_argument('--save', default='out', help='Prefix to save graphics')
parser.add_argument('--offsets', nargs='+', default=[], help='Vertical offsets for each molecule specified manually')
parser.add_argument('--padding', type=int, default=0, help='Horizontal padding to add between each molecule (in angstroms)')
parser.add_argument('--axes', action='store_true', default=False, help='Draw x and y axes')
parser.add_argument('--format', default='png', help='Format to save graphics', choices=['svg','pdf','png'])
parser.add_argument('--membrane', default=None, choices=[None, 'fancy', 'box'], help='Draw membrane on X axis')
parser.add_argument('--dpi', type=int, default=300, help='DPI to use if exporting to raster formats (i.e. PNG)')

# for simulating according to stoichiometry
parser.add_argument('--csv', help='Table of protein information')
parser.add_argument('--num_mol', type=int, help='Total number of molecules in the scene')
parser.add_argument('--order_by', default='input', choices=['input', 'random', 'height'], help='How to order proteins in scene')
parser.add_argument('--background', action='store_true', default=False, help='Add background plane using same frequencies')

args = parser.parse_args()

def draw_membrane(width, height=40):
    axs = plt.gca()
    membrane_box = mpatches.FancyBboxPatch(
        [-100, 0], 1.5*width, -1*height,
        boxstyle=mpatches.BoxStyle("Round", pad=0.02),
        facecolor='#DDD5C7', ec='#DDD5C7', alpha=1, zorder=2)
    axs.add_patch(membrane_box)

# cartoon of lipid bilayer
def draw_membrane_fancy(width, height=40, jitter=False):
    axs = plt.gca()
    head_radius = 4
    num_lipids = int(2*width/head_radius)

    # color scheme
    lipid_head_fc='#D6D1EF'
    lipid_tail_fc='#A3DCEF'
    membrane_box_fc='#C4E7EF'

    # membrane box background
    membrane_box = mpatches.FancyBboxPatch(
        [-100, -1*head_radius+10], 1.5*width, -1*height+head_radius,
        boxstyle=mpatches.BoxStyle("Round", pad=0.02),
        facecolor=membrane_box_fc, alpha=1, zorder=1)
    axs.add_patch(membrane_box)

    # draw each lipid pair
    for i in range(num_lipids):
        top_jitter_y = 10
        bot_jitter_y = 10
        if jitter:
            top_jitter_y += (np.random.random()-0.5)*2*head_radius
            bot_jitter_y += (np.random.random()-0.5)*2*head_radius
        axs.add_line(mlines.Line2D([i*head_radius*2-100, i*head_radius*2-100], [-4+top_jitter_y, -18+top_jitter_y], zorder=1.5, c=lipid_tail_fc, linewidth=head_radius*.7, alpha=1, solid_capstyle='round'))
        axs.add_line(mlines.Line2D([i*head_radius*2-100, i*head_radius*2-100], [-38+bot_jitter_y, -24+bot_jitter_y], zorder=1.5, c=lipid_tail_fc, linewidth=head_radius*.7, alpha=1, solid_capstyle='round'))
        axs.add_patch(mpatches.Circle((i*head_radius*2-100, -1*head_radius+top_jitter_y), head_radius, facecolor=lipid_head_fc, ec='k', linewidth=0.3, alpha=1, zorder=2))
        axs.add_patch(mpatches.Circle((i*head_radius*2-100, -1*height+bot_jitter_y), head_radius, facecolor=lipid_head_fc, ec='k', linewidth=0.3, alpha=1, zorder=2))

# list of protein polygons to draw
object_list = []
num_files = 0

if args.files:
    for path in args.files:
        with open(path,'rb') as f:
            data = pickle.load(f)
            object_list.append(data)
    num_files = len(args.files)
elif args.csv:
    protein_data = dict()
    with open(args.csv) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            (name, stoich, path) = (row['name'], int(row['stoichiometry']), row['file'])
            with open(path,'rb') as f:
                data = pickle.load(f)
                protein_data[name] = (stoich, data)
    num_files = len(protein_data)
else:
    sys.exit("No input files specified, see options with --help")

if len(args.offsets) > 0:
    assert(len(args.files) == len(args.offsets))
    y_offsets = list(map(float, args.offsets))
else:
    y_offsets = np.zeros(len(object_list))

if args.num_mol is None:
    num_mol = num_files
else:
    num_mol = args.num_mol

if args.csv:
    # total sum of protein counts
    protein_names = np.array(sorted(protein_data.keys()))
    protein_stoich = np.array([protein_data[p][0] for p in protein_names])
    sum_stoich = np.sum(protein_stoich)
    stoich_weights = protein_stoich / sum_stoich

    # protein copy number
    sampled_protein = np.random.choice(protein_names, size=num_mol, p=stoich_weights)
    object_list = [protein_data[p][1] for p in sampled_protein]

    # assemble objects for background
    if args.background:
        scaling_factor = 0.6
        sampled_protein = np.random.choice(protein_names, int(num_mol*1/scaling_factor), p=stoich_weights)
        background_object_list = [protein_data[p][1] for p in sampled_protein]

# sort proteins
if args.order_by == "random":
    np.random.shuffle(object_list)
elif args.order_by == "height":
    object_list = sorted(object_list, key=lambda x: x['height'], reverse=True)

# set font options
font_options = {'family':'Arial', 'weight':'normal', 'size':10}
matplotlib.rc('font', **font_options)

# set up plot
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

# draw molecules
w=0
for i, o in enumerate(object_list):
    # infer if --residues and use lighter line width
    if len(o['polygons']) > 50:
        lw=0.1
    else:
        lw=0.5
    for p in o['polygons']:
        xy = p.get_xy()
        axs.fill(xy[:,0]+w, xy[:,1], fc=p.get_facecolor(), ec='k', linewidth=lw, zorder=3)
    w += o['width']+args.padding

if args.background:
    background_w=0
    for i, o in enumerate(background_object_list):
        for p in o['polygons']:
            xy = p.get_xy()
            axs.fill((xy[:,0]+background_w)*scaling_factor, (xy[:,1])*scaling_factor, fc=p.get_facecolor(), linewidth=0, zorder=1)
        background_w += (o['width']+args.padding)

if args.membrane is not None:
    plt.gca().set_ylim((-45,350))
    if args.membrane == 'box':
        draw_membrane(width=w)
    elif args.membrane == 'fancy':
        draw_membrane_fancy(width=w)

plt.savefig(args.save+'.'+args.format, transparent=True, pad_inches=0, bbox_inches='tight', dpi=args.dpi)
