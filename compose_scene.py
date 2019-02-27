import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import lines, text, cm
import shapely.geometry as sg
import shapely.ops as so
import os, sys, re, argparse, csv, pickle
from scipy.signal import savgol_filter
from scipy import interpolate
import glob

parser = argparse.ArgumentParser(description='Structure SVG compositing',  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--files', nargs='+', help='Pickled objects to load')
parser.add_argument('--offsets', nargs='+', default=[], help='Vertical offsets for each molecule specified manually')
parser.add_argument('--axes', action='store_true', default=False, help='Draw x and y axes')
parser.add_argument('--format', default='png', help='Format to save graphics', choices=['svg','pdf','png'])
parser.add_argument('--membrane', action='store_true', default=False, help='Pickled objects to load')
args = parser.parse_args()

def draw_membrane(width):
    axs = plt.gca()
    membrane_box = mpatches.FancyBboxPatch(
        [-100, 0], 2*width, -40,
        boxstyle=mpatches.BoxStyle("Round", pad=0.02),
        edgecolor='none',facecolor='#DDD5C7',alpha=0.5, zorder=1)
    axs.add_patch(membrane_box)

# read in saved Python objects
object_list = []
for path in args.files:
    with open(path,'rb') as f:
        data = pickle.load(f)
        object_list.append(data)

if len(args.offsets) > 0:
    assert(len(args.files) == len(args.offsets))
    y_offsets = list(map(float, args.offsets))
else:
    y_offsets = np.zeros(len(object_list))

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
spacer=30
w=object_list[0]['width']/2
for i, o in enumerate(object_list):
    w += o['width']/2
    for p in o['polygons']:
        xy = p.get_xy()
        axs.fill(xy[:,0]+w, xy[:,1]+y_offsets[i], fc=p.get_facecolor(), ec='k', linewidth=0.5, zorder=2)
    w += o['width']/2 + spacer

if args.membrane:
    draw_membrane(width=w)

plt.savefig('compose.'+args.format, dpi=300)
