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

def rotation_matrix_2d(theta):
    return np.array([[np.cos(theta), -1*np.sin(theta)],[np.sin(theta), np.cos(theta)]])

def draw_object(o, axs, offset=[0,0], flip=False, background=False, scaling=1, zorder=3, recenter=None, color=None, linewidth=None):
    """
    add object to axes
    abstracting some of the drawing details, might want to move to a class
    color option overrides original color
    """

    # use thinner line with many objects (e.g. with --residues)
    if linewidth is None:
        if len(o['polygons']) > 50:
            lw=0.1
        else:
            lw=0.2
    else:
        lw = linewidth

    for p in o['polygons']:

        # get polygon coordinates and transform if necessary
        xy = p.get_xy()
        if recenter is not None:
            # optionally shift coordinates before rotation
            xy -= recenter
        if flip:
            # TODO easier to use object height/width if accurate
            xy = np.dot(xy, np.array([[-1,0],[0,-1]]))
            offset_x = np.min(xy[:,0])
            offset_y = np.min(xy[:,1])
            xy -= np.array([offset_x, offset_y])

        # either use existing color or recolor with global color scheme
        if color is not None:
            fc = color
        else:
            fc = p.get_facecolor()

        # fill polygon, preset for semi-transparent background layer
        if not background:
            axs.fill((xy[:,0]+offset[0])*scaling, (xy[:,1]+offset[1])*scaling, fc=fc, ec='k', linewidth=lw*scaling, zorder=zorder)
        else:
            axs.fill((xy[:,0]+offset[0])*scaling, (xy[:,1]+offset[1])*scaling, fc=fc, ec='k', alpha=0.8, linewidth=lw*scaling, zorder=1)

class membrane_cartoon:
    def __init__(self, width, thickness, axes, base_y=0):
        self.width = width
        self.thickness = thickness
        self.y = base_y
        self.axes = axes
        # other constants
        self.head_radius = 4

    def flat(self):
        self.height_at = lambda x: self.y + self.thickness/2

    def sinusoidal(self, frequency=1, amplitude=1):
        self.height_at = lambda x: self.y + self.thickness/2*amplitude*np.sin(x*frequency*2*np.pi/self.width)

    def interpolate(self, x, y, kind='linear'):
        #self.height_at = interpolate.interp1d(x, y, kind=kind)
        self.height_fn = interpolate.PchipInterpolator(x, y)
        self.height_at = lambda x: self.height_fn(x) + self.y

    def draw(self, lipids=False):
        membrane_box_fc='#C4E7EF'
        membrane_x = np.linspace(0,self.width,200)
        membrane_y_top = np.array([self.height_at(x) for x in membrane_x])
        membrane_y_bot = membrane_y_top-self.thickness
        plt.fill_between(membrane_x, membrane_y_top-self.head_radius, membrane_y_bot+self.head_radius, color=membrane_box_fc, zorder=1.6)
        if lipids:
            # color scheme
            lipid_head_fc='#D6D1EF'
            lipid_tail_fc='#A3DCEF'
            num_lipids = int(self.width/(2*self.head_radius))
            for i in range(num_lipids):
                membrane_y = self.height_at(i/num_lipids*self.width)
                self.axes.add_line(mlines.Line2D([i*self.head_radius*2, i*self.head_radius*2], [-4+membrane_y, -18+membrane_y], zorder=1.7, c=lipid_tail_fc, linewidth=self.head_radius*.7, alpha=1, solid_capstyle='round'))
                self.axes.add_line(mlines.Line2D([i*self.head_radius*2, i*self.head_radius*2], [-38+membrane_y, -24+membrane_y], zorder=1.7, c=lipid_tail_fc, linewidth=self.head_radius*.7, alpha=1, solid_capstyle='round'))
                self.axes.add_patch(mpatches.Circle((i*self.head_radius*2, -1*self.head_radius+membrane_y), self.head_radius, facecolor=lipid_head_fc, ec='k', linewidth=0.3, alpha=1, zorder=2))
                self.axes.add_patch(mpatches.Circle((i*self.head_radius*2, -1*self.thickness+membrane_y), self.head_radius, facecolor=lipid_head_fc, ec='k', linewidth=0.3, alpha=1, zorder=2))

def make_scene(args):
    # list of protein polygons to draw
    object_list = []
    num_files = 0

    if args.files:
        for path in args.files:
            with open(path,'rb') as f:
                data = pickle.load(f)
                object_list.append(data)
                coords = np.concatenate([p.get_xy() for p in data['polygons']])
                data['height'] = np.max(coords[:,1]) # temporary bugfix
        num_files = len(args.files)
    elif args.csv:
        protein_data = dict()
        with open(args.csv) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                (name, stoich, path) = (row['name'], float(row[args.sample_from]), row['file'])
                with open(path+'.pickle','rb') as f:
                    data = pickle.load(f)
                    data['name'] = name
                    data['stoichiometry'] = stoich
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
            scaling_factor = 0.7
            sampled_protein = np.random.choice(protein_names, int(num_mol*1/scaling_factor), p=stoich_weights)
            background_object_list = [protein_data[p][1] for p in sampled_protein]
    else:
        if 'name' in object_list[0]:
            protein_names = [o['name'] for o in object_list]
        else:
            for i,o in enumerate(object_list):
                o['name'] = i
            protein_names = range(len(object_list))

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

    if args.recolor:
        # default cmap is hsv. for discrete could try Set1 or Pastel1
        if len(args.recolor_cmap) == 1:
            cmap = cm.get_cmap(args.recolor_cmap[0])
        else:
            cmap = ListedColormap(args.recolor_cmap)
        color_scheme = dict()
        for i,c in enumerate(object_list):
            name = c['name']
            if isinstance(cmap, ListedColormap):
                color_scheme[name] = cmap(i)
            else:
                color_scheme[name] = cmap(i/len(object_list))

    if args.membrane is not None:
        total_width = np.sum([o['width'] for o in object_list])+len(object_list)*args.padding
        membrane = membrane_cartoon(width=total_width, axes=axs, thickness=40)

        # <-------------------------------------------------------------------------
        # testing new code
        if args.membrane_interface and len(object_list) > 2:
            skyline_pos = []
            w=0
            for i, o in enumerate(object_list):
                skyline_pos.append([w,o['height']+np.random.rand()*10])
                skyline_pos.append([w+o['width'],o['height']+np.random.rand()*10])
                w += o['width'] + args.padding
            skyline_pos = np.array(skyline_pos).T
            membrane2 = membrane_cartoon(width=total_width, axes=axs, thickness=40, base_y=50)
            membrane2.interpolate(skyline_pos[0], skyline_pos[1])
            membrane2.draw(lipids=True)
        # <-------------------------------------------------------------------------

        if args.membrane == "flat":
            membrane.flat()
        elif args.membrane == "arc":
            membrane.sinusoidal(frequency=0.5, amplitude=2)
        elif args.membrane == "wave":
            membrane.sinusoidal(frequency=2, amplitude=2)
        membrane.draw(lipids=args.membrane_lipids)

    # draw molecules
    w=0
    for i, o in enumerate(object_list):
        y_offset = membrane.height_at(w+o['bottom'][0])-10
        if args.recolor:
            color = color_scheme[o['name']]
        else:
            color = None
        draw_object(o, axs, offset=[w, y_offset], color=color)
        w += o['width']+args.padding

    if args.background:
        background_w=0
        for i, o in enumerate(background_object_list):
            draw_object(o, axs, offset=[background_w, 0], scaling=scaling_factor, background=True)
            background_w += (o['width']+args.padding)

    plt.savefig(args.save+'.'+args.format, transparent=True, pad_inches=0, bbox_inches='tight', dpi=args.dpi)