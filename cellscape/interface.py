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

from .scene import draw_object

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
        self.axes.fill_between(membrane_x, membrane_bot_y, membrane_bot_y-self.thickness, color=color, zorder=1.6, capstyle='round', joinstyle='miter')

        # plot top membrane
        self.axes.fill_between(membrane_x, membrane_top_y, membrane_top_y+self.thickness, color=color, zorder=1.6, capstyle='round', joinstyle='round')

def plot_pairs(pairs, thickness=40, padding=50, align="bottom", membrane_color="#E8E8E8", colors=None, axes=True, linewidth=None):
    fig, axs = plt.subplots()
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

    assert(align in ["top","bottom","middle"])

    # get all the interface heights
    all_heights = np.array([p[0]['height']+p[1]['height'] for p in pairs])
    max_height = np.max(all_heights)

    # calculate membrane geometry
    bot_y = []
    top_y = []
    lengths = []
    for p in pairs:
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

    # draw membrane
    mem = MembraneInterface(axes=axs, lengths=lengths, bottom_y=bot_y, top_y=top_y, padding=padding, thickness=thickness)
    mem.draw(color=membrane_color)

    # draw proteins
    w=0
    for i, o in enumerate(pairs):
        o1, o2 = o
        this_width = max(o1['width'], o2['width'])
        this_height = o1['height']+o2['height']
        y_offset = bot_y[i]
        if colors is not None:
            color_top = colors[0][i]
            color_bot = colors[1][i]
        else:
            color_top = None
            color_bot = None
        draw_object(o1, axs, offset=[w+(this_width-o1['width'])/2, y_offset], color=color_top, linewidth=linewidth)
        draw_object(o2, axs, offset=[w+(this_width-o2['width'])/2, y_offset+o1['height']], flip=True, color=color_bot, linewidth=linewidth)
        w += this_width+padding

    fig.set_size_inches(18.5, 10.5)
    return fig
