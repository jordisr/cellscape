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
parser.add_argument('--grid', action='store_true', help='Grid display')
args = parser.parse_args()

#object_filenames = sys.argv[1:]
object_list = []
for path in glob.glob("*.pickle"):
    with open(path,'rb') as f:
        data = pickle.load(f)
        object_list.append(data)

object_list = sorted(object_list, key=lambda x: x['height'])

cmap = cm.get_cmap('viridis')
cmap_x = np.linspace(0.0,1.0,len(object_list))
sequential_colors = [cmap(x) for x in cmap_x]

fig, axs = plt.subplots()
axs.set_aspect('equal')
plt.axis('off')

if args.grid:
    i_window = int(len(object_list)/2)
    for i in range(0,len(object_list),i_window):
        cum_width = 0
        for j,obj in enumerate(object_list[i:i+i_window]):
            obj_x = cum_width
            obj_y = i*30
            print(obj['name'], obj['width'])
            outline = obj['outline']
            outline[0] -= np.mean(outline[0])
            outline[1] -= np.mean(outline[1])
            axs.fill(outline[0]+np.array([obj_x]), outline[1]+np.array([obj_y]), alpha=1, fc=sequential_colors[i % len(sequential_colors)], ec='k',zorder=2,linewidth=0.5)
            cum_width += obj['width']*1.5

    fig = plt.gcf()
    fig.set_size_inches(30, 10.5)
    fig.savefig('scene.pdf',transparent=True, pad_inches=0, bbox_inches='tight')

else:
    membrane_box = mpatches.FancyBboxPatch(
        [-100, 0], len(object_list)*100, -40,
        boxstyle=mpatches.BoxStyle("Round", pad=0.02),
        edgecolor='none',facecolor='#DDD5C7',alpha=0.5)
    axs.add_patch(membrane_box)

    cum_width = 0
    for i,obj in enumerate(object_list):
        obj_x = cum_width
        obj_y = 0
        print(obj['name'], obj['width'])
        if 'domains' in obj:
            for i, path in enumerate(obj['domain_paths']):
                axs.fill(path[0]+np.array([obj_x]), path[1]+np.array([obj_y]), alpha=1, fc=sequential_colors[i % len(sequential_colors)], ec='k',zorder=2,linewidth=0)
        else:
            outline = obj['outline']
            #outline[0] -= np.mean(outline[0])
            #outline[1] -= np.mean(outline[1])
            axs.fill(outline[0]+np.array([obj_x]), outline[1]+np.array([obj_y]), alpha=1, fc=sequential_colors[i % len(sequential_colors)], ec='k',zorder=2,linewidth=0.5)
        cum_width += obj['width']*1.5

    fig = plt.gcf()
    fig.set_size_inches(30, 10.5)
    fig.savefig('scene.pdf',transparent=True, pad_inches=0, bbox_inches='tight')
