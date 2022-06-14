"""
Testing code for visualizing protein interactions across membrane interfaces
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import lines, text, cm
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from scipy import interpolate
import shapely.geometry as sg
import shapely.ops as so
import os, sys, argparse, pickle
import glob
import csv

from cellscape.cartoon import plot_polygon, shade_from_color

class MembraneInterface:
    """
    just use piecemeal flat + connector, no lipids
    """
    def __init__(self, axes, lengths, bottom_y, top_y, thickness=40, padding=10, base_y=0):

        # axes
        self.axes = axes

        # membrane thickness (angstroms)
        self.thickness = thickness

        # padding between each segment, scalar
        self.padding = padding

        # length of each segment, array
        self.lengths = lengths

        # y coordinate of each bottom membrane segment, array
        self.bottom_y = bottom_y

        # y coordinate of each top membrane segment, array
        self.top_y = top_y

        assert(len(lengths) == len(top_y))
        assert(len(top_y) == len(bottom_y))

    def draw(self, color='#C4E7EF'):
        if isinstance(color, (list,tuple)):
            top_color = color[0]
            bot_color = color[1]
        else:
            top_color = color
            bot_color = color

        membrane_x = []
        membrane_bot_y = []
        membrane_top_y = []
        x_cum = 0
        for i, w in enumerate(self.lengths):
            membrane_x.append(x_cum)
            x_cum += w
            membrane_x.append(x_cum)
            x_cum += self.padding

            membrane_bot_y.append(self.bottom_y[i])
            membrane_bot_y.append(self.bottom_y[i])

            membrane_top_y.append(self.top_y[i])
            membrane_top_y.append(self.top_y[i])

        membrane_x = np.array(membrane_x)
        membrane_bot_y = np.array(membrane_bot_y)
        membrane_top_y = np.array(membrane_top_y)

        # plot bottom membrane
        self.axes.fill_between(membrane_x, membrane_bot_y, membrane_bot_y-self.thickness, color=bot_color, zorder=1.6, capstyle='round', joinstyle='miter')

        # plot top membrane
        self.axes.fill_between(membrane_x, membrane_top_y, membrane_top_y+self.thickness, color=top_color, zorder=1.6, capstyle='round', joinstyle='round')

def plot_pairs(pairs, labels=None, thickness=40, padding=50, align="bottom", membrane_color="#E8E8E8", colors=None, axes=True, linewidth=None, sort=False):

    assert align in ["bottom", "middle", "top"]
    assert sort in [False, "height", "horseshoe"]

    if labels is not None:
        assert len(labels) == len(pairs)

    # optionally sort proteins by height
    if sort == "height":
        pair_heights = np.array(list(map(lambda x: x[0]['height']+x[1]['height'], pairs)))
        sorted_order = np.argsort(pair_heights)[::-1]
        pairs_ = [pairs[i] for i in sorted_order]
        if labels is not None:
            labels_ = [labels[i] for i in sorted_order]

    elif sort == "horseshoe":
        pair_heights = np.array(list(map(lambda x: x[0]['height']+x[1]['height'], pairs)))
        sorted_order = np.argsort(pair_heights)[::-1]
        new_order = np.zeros_like(sorted_order)
        first_half = sorted_order[::2]
        if len(sorted_order) % 2:
            second_half = sorted_order[-2::-2]
        else:
            second_half = sorted_order[::-2]
        new_order[:len(first_half)] = first_half
        new_order[len(first_half):] = second_half

        pairs_ = [pairs[i] for i in new_order]
        if labels is not None:
            labels_ = [labels[i] for i in new_order]
        if colors is not None:
            colors_ = [colors[i] for i in new_order]

    else:
        pairs_ = pairs[:]
        if labels is not None:
            labels_ = labels[:]

    fig, axs = plt.subplots(figsize=(11,8.5))
    axs.set_aspect('equal')

    if axes:
        axs.xaxis.grid(False)
        axs.yaxis.grid(False)
        axs.axes.xaxis.set_ticklabels([])
    else:
        plt.axis('off')
        plt.gca().set_axis_off()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
        plt.margins(0,0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())

    # set font options
    font_options = {'family':'Arial', 'weight':'normal', 'size':10}
    matplotlib.rc('font', **font_options)

    assert(align in ["top","bottom","middle"])

    # get all the interface heights
    all_heights = np.array([p[0]['height']+p[1]['height'] for p in pairs_])
    max_height = np.max(all_heights)

    # calculate membrane geometry
    bot_y = []
    top_y = []
    lengths = []
    for p in pairs_:
        o1, o2 = p
        if align == "bottom":
            top_y.append(o1['height']+o2['height'])
            bot_y.append(0)
        elif align == "top":
            top_y.append(max_height)
            bot_y.append(max_height-(o1['height']+o2['height']))
        elif align == "middle":
            top_y.append(max_height-(max_height-o1['height']-o2['height'])/2)
            bot_y.append((max_height-o1['height']-o2['height'])/2)
        lengths.append(max(o1['width'], o2['width']))
    top_y = np.array(top_y)
    bot_y = np.array(bot_y)
    lengths = np.array(lengths)

    total_width = np.sum(lengths)+len(pairs_)*padding

    # draw membrane
    mem = MembraneInterface(axes=axs, lengths=lengths, bottom_y=bot_y, top_y=top_y, padding=padding, thickness=thickness)
    mem.draw(color=membrane_color)

    # draw proteins
    w=0
    for i, o in enumerate(pairs_):
        o1, o2 = o
        this_width = max(o1['width'], o2['width'])
        this_height = o1['height']+o2['height']
        y_offset = bot_y[i]
        if colors is not None:
            color_top = colors_[i][1]
            color_bot = colors_[i][0]
        else:
            color_top = None
            color_bot = None

        # TODO rotation needs to be a little cleaner, making some assumptions here
        for p in o1["polygons"]:
            xy = np.array(o1['polygons'][0]['polygon'].exterior.xy) # assuming first polygon is outline, TODO fix
            recenter = np.array([np.min(xy[:,0]), np.min(xy[:,1])])
            facecolor = shade_from_color(color_bot, p.get("shade", 0.5), range=p.get("shading_range", 0.4))
            plot_polygon(p['polygon'], axes=axs, translate_pre=[-recenter[0]+w+(this_width-o1['width'])/2, -recenter[1]+y_offset], flip=False, facecolor=facecolor, linewidth=p['linewidth'])

        for p in o2["polygons"]:
            xy = np.array(o2['polygons'][0]['polygon'].exterior.xy) # assuming first polygon is outline, TODO fix
            recenter = np.array([np.min(xy[:,0]), np.min(xy[:,1])])
            facecolor = shade_from_color(color_top, p.get("shade", 0.5), range=p.get("shading_range", 0.4))
            plot_polygon(p['polygon'], axes=axs, translate_pre=-1*recenter, translate_post=[w+(this_width+o2['width'])/2, y_offset+10+o2['height']+o1['height']], flip=True, facecolor=facecolor, linewidth=p['linewidth'])
        #plot_polygon(o2["polygons"][0]['polygon'], axes=axs, offset=[w+(this_width-o2['width'])/2, y_offset+o1['height']], flip=True, facecolor=o2["polygons"][0]['facecolor'], linewidth=linewidth)

        if labels_ is not None:
            angstroms_per_inch = total_width/11
            fontsize = total_width*0.5/len(pairs_)/angstroms_per_inch*72
            font_inches = fontsize/72
            plt.text(w+this_width/2,  y_offset+this_height+50, labels_[i][1], rotation=90, fontsize=fontsize, va='bottom', ha='center')
            plt.text(w+this_width/2, y_offset-1.1*angstroms_per_inch*font_inches, labels_[i][0], rotation=90, fontsize=fontsize, va='top', ha='center')

        w += this_width+padding

    fig.set_size_inches(18.5, 10.5)
    return fig
