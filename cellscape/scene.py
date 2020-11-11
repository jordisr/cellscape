import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import lines, text, cm
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from scipy import interpolate
import os, sys, argparse, pickle
import glob
import csv

from .cartoon import plot_polygon, shade_from_color, placeholder_polygon

def rotation_matrix_2d(theta):
    """Return matrix to rotate 2D coordinates by angle theta."""
    return np.array([[np.cos(theta), -1*np.sin(theta)],[np.sin(theta), np.cos(theta)]])

def draw_object(o, axs, offset=[0,0], flip=False, background=False, scaling=1, zorder=3, recenter=None, color=None, linewidth=None):
    # older function for adding proteins to membrane visualization
    # now just use plot_polygon

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

class Membrane:
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

        membrane_x = np.linspace(0,self.width,200)
        membrane_y_top = np.array([self.height_at(x) for x in membrane_x])
        membrane_y_bot = membrane_y_top-self.thickness

        if lipids:
            membrane_box_fc='#C4E7EF'
            lipid_head_fc='#D6D1EF'
            lipid_tail_fc='#A3DCEF'
            plt.fill_between(membrane_x, membrane_y_top-self.head_radius, membrane_y_bot+self.head_radius, color=membrane_box_fc, zorder=1.6)
            num_lipids = int(self.width/(2*self.head_radius))
            for i in range(num_lipids):
                membrane_y = self.height_at(i/num_lipids*self.width)
                self.axes.add_line(mlines.Line2D([i*self.head_radius*2, i*self.head_radius*2], [-4+membrane_y, -18+membrane_y], zorder=1.7, c=lipid_tail_fc, linewidth=self.head_radius*.7, alpha=1, solid_capstyle='round'))
                self.axes.add_line(mlines.Line2D([i*self.head_radius*2, i*self.head_radius*2], [-38+membrane_y, -24+membrane_y], zorder=1.7, c=lipid_tail_fc, linewidth=self.head_radius*.7, alpha=1, solid_capstyle='round'))
                self.axes.add_patch(mpatches.Circle((i*self.head_radius*2, -1*self.head_radius+membrane_y), self.head_radius, facecolor=lipid_head_fc, ec='k', linewidth=0.3, alpha=1, zorder=2))
                self.axes.add_patch(mpatches.Circle((i*self.head_radius*2, -1*self.thickness+membrane_y), self.head_radius, facecolor=lipid_head_fc, ec='k', linewidth=0.3, alpha=1, zorder=2))

        else:
            membrane_box_fc='#C8C8C8'
            plt.fill_between(membrane_x, membrane_y_top, membrane_y_bot, color=membrane_box_fc, zorder=1.6)

def make_scene(args):
    """Build a scene in one-go. Called when running ``cellscape scene``."""

    assert args.save.split('.')[-1] in ['png','pdf','svg','ps'], "image format not recognized"

    # list of protein polygons to draw
    object_list = []
    num_files = 0

    # set random seed for reproducibility
    if args.seed:
        np.random.seed(args.seed)

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
                (name, stoich, path) = (row['name'], float(row[args.sample_from]), row.get('file'))
                if path is not "":
                    with open(path,'rb') as f:
                        data = pickle.load(f)
                        data['name'] = name
                        data['stoichiometry'] = stoich
                        # TEST specifying color in CSV file
                        if 'color' in row:
                            data['color'] = row['color']
                        protein_data[name] = (stoich, data)
                elif args.use_placeholders:
                    height = float(row.get('height'))*10 # assuming in nanometers
                    data = {'name':name, 'stoichiometry':stoich, 'height':height, 'bottom':np.array([25,0]), 'width':50, 'polygons':[{'polygon':placeholder_polygon(height), 'edgecolor':'k', 'linewidth':1, 'facecolor':"#eeeeee"}]}
                    protein_data[name] = (stoich, data)

        num_files = len(protein_data)

    else:
        sys.exit("No input files specified, see options with --help")

    if len(args.offsets) > 0:
        assert(len(args.files) == len(args.offsets))
        y_offsets = list(map(float, args.offsets))
    else:
        y_offsets = np.zeros(len(object_list))

    if args.csv:
        # total sum of protein counts
        protein_names = np.array(list(protein_data.keys()))
        protein_stoich = np.array([protein_data[p][0] for p in protein_names])
        sum_stoich = np.sum(protein_stoich)
        stoich_weights = protein_stoich / sum_stoich

        if args.num_mol > 0:
            # protein copy number
            sampled_protein = np.random.choice(protein_names, size=args.num_mol, p=stoich_weights)
            object_list = [protein_data[p][1] for p in sampled_protein]
        else:
            object_list = [protein_data[p][1] for p in protein_names]

        # assemble objects for background
        if args.background and args.num_mol > 0:
            scaling_factor = 0.7
            sampled_protein = np.random.choice(protein_names, int(args.num_mol*1/scaling_factor), p=stoich_weights)
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
    elif args.order_by == "top":
        # TODO should be renamed, maybe length for overall size and height for above membrane?
        object_list = sorted(object_list, key=lambda x: x['top'][1], reverse=True)

    # set font options
    font_options = {'family':'Arial', 'weight':'normal', 'size':10}
    matplotlib.rc('font', **font_options)

    # set up plot
    # POSSIBLE BUG: while coordinates and scale are prserved in the pickle files,
    # this doesn't necessarily apply to the images. Hence if someone tries to
    # manually add a protein to a scene that has been generated there could be
    # sizing issues. Is the solution to use  a constant angstrom/inch scaling?
    scene_height_in = args.fig_height
    scene_width_in = args.fig_width
    fig, axs = plt.subplots(figsize=(scene_width_in, scene_height_in))
    axs.set_aspect('equal')

    if args.axes:
        plt.axis('on')
        axs.xaxis.grid(False)
        axs.yaxis.grid(True)
        axs.axes.xaxis.set_ticklabels([])
        axs.autoscale()
        plt.margins(0.01,0.01) # is this needed?

    else:
        plt.axis('off')
        plt.gca().set_axis_off()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
        axs.autoscale()
        plt.margins(0.01,0.01)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())

    if args.recolor:
        # default cmap is hsv. for discrete could try Set1 or Pastel1
        if len(args.recolor_cmap) == 1:
            cmap = cm.get_cmap(args.recolor_cmap[0])
        else:
            # TESTING interpret as continous color scheme
            # cmap = LinearSegmentedColormap.from_list("cmap", args.recolor_cmap)
            cmap = ListedColormap(args.recolor_cmap)
        color_scheme = dict()
        for i,c in enumerate(sorted(object_list, key=lambda x: x['height'])):
            name = c['name']
            if isinstance(cmap, ListedColormap):
                color_scheme[name] = cmap(i)
            else:
                color_scheme[name] = cmap(i/len(object_list))

    # TESTING
    # so colors are by height (what about duplicated molecules)
    # np.random.shuffle(object_list)

    total_width = np.sum([o['width'] for o in object_list])+len(object_list)*args.padding
    if args.membrane is not None:
        membrane = Membrane(width=total_width, axes=axs, thickness=args.membrane_thickness)

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
        if args.membrane is not None and not args.no_membrane_offset:
            y_offset = membrane.height_at(w+o['bottom'][0])-10
        else:
            y_offset = 0
        for p in o["polygons"]:
            # TODO change to dict.get() call to have default
            if args.recolor:
                if 'color' in o:
                    # check if color specified in CSV file
                    facecolor = o['color']
                    edgecolor = p["edgecolor"]
                else:
                    # use color scheme from recolor_cmap
                    facecolor = color_scheme[o['name']]
                    edgecolor = 'black'
                if "shade" in p:
                    # TODO export shading_range from polygons as well
                    facecolor = shade_from_color(facecolor, p["shade"], range=p.get("shading_range", 0.4)) # using default from cartoon.py, could change
            else:
                # use color already specified
                facecolor = p["facecolor"]
                edgecolor = p["edgecolor"]

            plot_polygon(p["polygon"], translate_pre=[w, y_offset], facecolor=facecolor, edgecolor=edgecolor, linewidth=p["linewidth"], zorder_mod=p.get("zorder", 0))
            if args.labels:
                # option is experimental, text needs to be properly sized and placed
                # testing use of figure width in inches (specified above) and total width in angstroms to infer appropriate font size
                #plt.text(w+o['width']/2,-100, o.get("name", ""), rotation=90, fontsize=fontsize)
                # 1.1 and 0.6 numbers chosen through experimentation, best way would be to look at length of labels in characters
                angstroms_per_inch = total_width/scene_width_in
                fontsize = total_width*args.label_size/len(object_list)/angstroms_per_inch*72
                font_inches = fontsize/72
                # TODO better text positioning, allow for top/bottom selection
                if args.label_orientation == "vertical":
                    plt.text(w+o['width']/2,o['bottom'][1]-1.1*angstroms_per_inch*font_inches, o.get("name", ""), rotation=90, fontsize=fontsize, va='top', ha='center') # vertical text (below)
                elif args.label_orientation == "horizontal":
                    plt.text(w+o['width']/2,o['top'][1]+2*angstroms_per_inch*font_inches, o.get("name", ""), rotation=0, fontsize=fontsize, va='top', ha='center') # horizontal text (above)
                elif args.label_orientation == "diagonal":
                    plt.text(w+o['width']/5,o['top'][1]+angstroms_per_inch*font_inches, o.get("name", ""), rotation=45, fontsize=fontsize) # diagonal text (above)
        w += o['width']+args.padding

    if args.background:
        background_w=0
        for i, o in enumerate(background_object_list):
            # draw_object(o, axs, offset=[background_w, 0], scaling=scaling_factor, background=True)
            for p in o["polygons"]:
                plot_polygon(p["polygon"], offset=[background_w, 0], scale=scaling_factor, zorder_mod=p.get("zorder", -2), facecolor=p["facecolor"], edgecolor=p["edgecolor"], linewidth=p["linewidth"]*scaling_factor)
            background_w += (o['width']+args.padding)

    plt.savefig(args.save, transparent=True, pad_inches=0, bbox_inches='tight', dpi=args.dpi)
