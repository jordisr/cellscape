import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import shapely.geometry as sg
import shapely.ops as so
import pickle
import os
import sys
import colorsys
from Bio.PDB import *

from .util import *

def safe_union_accumulate(polys):
    # union of list of shapely polygons
    # maybe slower than cascaded_union/unary_union but catches topology errors
    # TODO figure out a better way of handling these exceptions
    exception_counter = 0
    u = polys[0]
    for p in polys:
        try:
            u = u.union(p)
        except:
            exception_counter += 1
            pass
    return u

def safe_union(a, b):
    # catch topology errors
    try:
        c = a.union(b)
    except:
        return a
    return c

def scale_line_width(x, lw_min, lw_max):
    return lw_max*(1-x) + lw_min*x

def shade_from_color(color, x, range):
    (r, g, b, a) = mcolors.to_rgba(color)
    h, l, s = colorsys.rgb_to_hls(r,g,b)
    l_dark = max(l-range/2, 0)
    l_light = min(l+range/2, 1)
    l_new = l_dark*(1-x) + l_light*x
    return colorsys.hls_to_rgb(h, l_new, s)

def get_sequential_colors(colors='Set1', n=1):
    """
    Sample n colors sequentially from a named matplotlib ColorMap.
    """
    # uses matplotlib.colors.ColorMap.N to distinguish continuous/discrete
    cmap = matplotlib.cm.get_cmap(colors)
    if cmap.N == 256:
        # continuous color map
        sequential_colors = [cmap(x) for x in np.linspace(0.0,1.0, n)]
    else:
        # discrete color map
        sequential_colors = [cmap(x) for x in range(n)]
    return sequential_colors

def smooth_polygon(p, level=0):
    # somewhat arbitrary but a lot easier than interpolation
    if level == 0:
        return p.simplify(0.3).buffer(-2, join_style=1).buffer(3, join_style=1)
    elif level == 1:
        return p.simplify(1).buffer(3, join_style=1).buffer(-5, join_style=1).buffer(4, join_style=1)
    elif level == 2:
        return p.simplify(3).buffer(5, join_style=1).buffer(-9, join_style=1).buffer(5, join_style=1)
    elif level == 3:
        return p.simplify(0.1).buffer(2, join_style=1)
    else:
        return p

def ring_coding(ob):
    # https://sgillies.net/2010/04/06/painting-punctured-polygons-with-matplotlib.html
    # The codes will be all "LINETO" commands, except for "MOVETO"s at the
    # beginning of each subpath
    #n = len(ob.coords)
    n = len(np.asarray(ob))
    codes = np.ones(n, dtype=Path.code_type) * Path.LINETO
    codes[0] = Path.MOVETO
    return codes

def placeholder_polygon(height, buffer_width=25, origin=[0,0]):
    return sg.LineString([(buffer_width+origin[0],0+origin[1]),(buffer_width+origin[0],height+origin[1])]).buffer(buffer_width)

def composite_polygon(cartoon, height_before, height_after, buffer_width=25):
    # placeholder + structure cartoon + placeholder
    if height_before > 0:
        before_poly =  placeholder_polygon(height_before, origin=cartoon.bottom_coord[:2]-[buffer_width,height_before], buffer_width=buffer_width)
        cartoon._styled_polygons.append({"polygon":before_poly, "facecolor":"#eeeeee", "shade":0.5, "edgecolor":'black', "linewidth":1, "zorder":-1})

    if height_after > 0:
        after_poly =  placeholder_polygon(height_after, origin=cartoon.top_coord[:2]-[buffer_width, 0], buffer_width=buffer_width)
        cartoon._styled_polygons.append({"polygon":after_poly, "facecolor":"#eeeeee", "shade":0.5, "edgecolor":'black', "linewidth":1, "zorder":-1})

    cartoon.image_height = cartoon.image_height + buffer_width + height_before + height_after
    cartoon.bottom_coord = cartoon.bottom_coord - np.array([0,height_before,0])
    cartoon.top_coord = cartoon.top_coord + np.array([0,height_after,0])

def export_placeholder(height, name, fname, buffer_width=25):
    # placeholder by itself
    poly =  placeholder_polygon(height, origin=[buffer_width, 0], buffer_width=buffer_width)
    styled_polygons = [{"polygon":poly, "facecolor":"#eeeeee", "shade":0.5, "edgecolor":'black', "linewidth":1, "zorder":-1}]

    data = {'polygons':styled_polygons, 'name':name, 'width':buffer_width*2, 'height':height+buffer_width, 'start':np.array([buffer_width,0]), 'end':np.array([height+2*buffer_width,0]), 'bottom':np.array([buffer_width,0]), 'top':np.array([height+2*buffer_width,0])}

    with open('{}.pickle'.format(fname),'wb') as f:
        pickle.dump(data, f)

def transform_coord(xy, translate_post=np.array([0,0]), translate_pre=np.array([0,0]), scale=1.0, flip=False):
    # 2d coordinates
    xy_ = xy
    if translate_pre is not None:
        # optionally shift coordinates before rotation
        xy_ += translate_pre
    if flip:
        xy_ = np.dot(xy_, np.array([[-1,0],[0,-1]]).T)
        #xy_ = np.dot(xy_, np.array([[-1,0],[0,-1]]))
        #offset_x = np.min(xy_[:,0])
        #offset_y = np.min(xy_[:,1])
        #xy_ -= np.array([offset_x, offset_y])
    return (xy_+translate_post)*scale

def polygon_to_path(polygon, min_interior_length=40, translate_pre=np.array([0,0]), translate_post=np.array([0,0]), scale=1.0, flip=False):
    # generate matplotlib Path object from Shapely polygon
    # filter out small interior holes and apply a scaling factor if desired
    #
    # https://sgillies.net/2010/04/06/painting-punctured-polygons-with-matplotlib.html
    # Convert coordinates to path vertices. Objects produced by Shapely's
    # analytic methods have the proper coordinate order, no need to sort.
    interiors = list(filter(lambda x: x.length > min_interior_length, polygon.interiors))
    vertices = np.concatenate(
                    [np.asarray(polygon.exterior)]
                    + [np.asarray(r) for r in interiors])
    codes = np.concatenate(
                [ring_coding(polygon.exterior)]
                + [ring_coding(r) for r in interiors])
    transformed_vertices = transform_coord(vertices, translate_pre=translate_pre, translate_post=translate_post, scale=scale, flip=flip)
    return Path(transformed_vertices, codes)

def plot_polygon(poly, facecolor='orange', edgecolor='k', linewidth=0.7, axes=None, zorder_mod=0, translate_pre=np.array([0,0]), translate_post=np.array([0,0]), scale=1.0, flip=False, recenter=None, min_area=7, linestyle='solid'):
    """Draw a Shapely polygon using matplotlib Patches."""
    if axes is None:
        axs = plt.gca()
        axs.set_aspect('equal')
    else:
        axs = axes
    if isinstance(poly, sg.polygon.Polygon):
        if poly.area > min_area:
            path = polygon_to_path(poly, translate_pre=translate_pre, translate_post=translate_post, scale=scale, flip=flip)
            patch = PathPatch(path, facecolor=facecolor, edgecolor='black', linewidth=linewidth, zorder=3+zorder_mod, linestyle=linestyle)
            axs.add_patch(patch)
    elif isinstance(poly, sg.multipolygon.MultiPolygon):
        for p in poly:
            plot_polygon(p, axes=axs, facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, scale=scale, zorder_mod=zorder_mod, translate_pre=translate_pre, translate_post=translate_post, flip=flip)

class Cartoon:
    """"""
    def __init__(self, polygons, residues, outline_by, back_outline, group_outlines, num_groups, dimensions, groups):
        # TODO currently just copying over all variables needed, should condense a little
        self._polygons = polygons
        self.residues_flat = residues
        self.outline_by = outline_by
        self.num_groups = num_groups
        self._back_outline = back_outline
        self._group_outlines = group_outlines
        self.dimensions = dimensions
        self.groups = groups

    def plot(self, colors=None, axes_labels=False, color_residues_by=None, edge_color="black", line_width=0.7,
        depth_shading=False, depth_lines=False, shading_range=0.4, smoothing=False, do_show=True, axes=None, save=None, dpi=300, placeholder=None):
        """
        Style and draw protein cartoon generated from the ``outline`` function.

        Mirroring biopython's phylogeny drawing options
        https://biopython.org/DIST/docs/api/Bio.Phylo._utils-module.html
        can optionally pass a matplotlib Axes instance instead of creating a new one
        if do_show is false then return axes object

        colors -- color scheme for plotting acceptable arguments are
            - named matplotlib-compatible color e.g. "red" (string)
            - hexadecimal color e.g. "#F8F8FF" (string)
            - list/tuple of colors e.g. ["red", "#F8F8FF"] (list/tuple)
            - dict of names to colors e.g. {"domain A": "red", "domain B":"blue"} (dict)
            - named discrete or continuous color scheme e.g. "Set1" (string)
        """

        self._styled_polygons = []

        if axes is None:
            # create a new matplotlib figure if none provided
            fig, axs = plt.subplots()
        else:
            assert(isinstance(axes, matplotlib.axes.Axes))
            axs = axes
    
        if axes_labels:
            axs.axis('on')
            axs.set_axis_on()
            axs.xaxis.grid(False)
            axs.yaxis.grid(True)
            axs.axes.xaxis.set_ticklabels([])
        else:
            axs.axis('off')
            axs.set_axis_off()
            #plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
            #axs.xaxis.set_major_locator(plt.NullLocator())
            #axs.yaxis.set_major_locator(plt.NullLocator())
    
        axs.set_aspect('equal')
        axs.margins(0,0)
        self._axes= axs

        # color schemes
        default_color = 'tab:blue'
        default_cmap = 'Set1'
        named_colors = [*mcolors.BASE_COLORS.keys(), *mcolors.TABLEAU_COLORS.keys(), *mcolors.CSS4_COLORS.keys(), *mcolors.XKCD_COLORS.keys()]

        # if outlining residues don't know number of color groups until plot is called
        if self.outline_by == "residue":
            if color_residues_by is None:
                num_colors_needed = 1
                residue_color_groups = {"all":self.residues_flat}
            else:
                residue_color_groups = group_by(self.residues_flat, lambda x: x.get(color_residues_by))
                num_colors_needed = len(residue_color_groups)
            self.num_groups = num_colors_needed

        # parse options and get list of base colors needed for plotting
        if colors is None:
            # choose default sequential color scheme based on number of colors needed
            if self.num_groups == 1:
                sequential_colors = [default_color]
            elif self.num_groups <= 9:
                sequential_colors = get_sequential_colors(colors="Set1", n=self.num_groups)
            elif self.num_groups <= 10:
                sequential_colors = get_sequential_colors(colors="tab10", n=self.num_groups)
            else:
                sequential_colors = get_sequential_colors(colors="tab20", n=self.num_groups)
        else:
            if isinstance(colors, dict):
                sequential_colors = []
            else:
                if isinstance(colors, str):
                    if self.num_groups == 1:
                        sequential_colors = [colors]
                    else:
                        sequential_colors = get_sequential_colors(colors=colors, n=self.num_groups)
                elif isinstance(colors, (list, tuple)):
                    if self.num_groups == 1:
                        if (len(colors) == 4) or (len(colors) == 3):
                            # assume single RGBA or RGB color
                            sequential_colors = [colors]
                        else:
                            sequential_colors = [colors[0]]
                    elif self.num_groups == len(colors):
                        sequential_colors = colors
                    else:
                        sys.exit("Insufficient colors provided")
        assert(len(sequential_colors) == self.num_groups)

        # color scheme represented as dict that maps group names to colors
        if self.outline_by == "residue":
            if len(sequential_colors) > 0:
                color_map = {k:sequential_colors[i] for i,k in enumerate(residue_color_groups.keys())}
            else:
                color_map = colors
        elif self.outline_by == "all":
            color_map = {None:sequential_colors[0]}
        else:
            if len(sequential_colors) > 0:
                color_map = {k:sequential_colors[i] for i,k in enumerate(self.groups)}
            else:
                color_map = colors
        assert(isinstance(color_map, dict))

        if self._back_outline is not None:
            if smoothing:
                smoothed_poly = smooth_polygon(self._back_outline, level=3)
                plot_polygon(smoothed_poly, facecolor="None", scale=1.0, axes=axs, edgecolor=edge_color, linewidth=2, zorder_mod=-1)
                self._styled_polygons.append({"polygon":smoothed_poly, "facecolor":"None", "edgecolor":edge_color, "linewidth":2})
            else:
                plot_polygon(self._back_outline, facecolor="None", scale=1.0, axes=axs, edgecolor=edge_color, linewidth=2, zorder_mod=-1)
                self._styled_polygons.append({"polygon":self._back_outline, "facecolor":"None", "edgecolor":edge_color, "linewidth":2})

        if len(self._group_outlines) > 0:
            for p in self._group_outlines:
                plot_polygon(p, facecolor="None", scale=1.0, axes=axs, edgecolor=edge_color, linewidth=1, zorder_mod=2)
                self._styled_polygons.append({"polygon":p, "facecolor":"None", "edgecolor":edge_color, "linewidth":1, "zorder":2})

        # TODO optionally show placeholder for unstructured regions
        if placeholder is not None:
            placeholder_poly = placeholder_polygon(placeholder-self.image_height, origin=[self.image_width/2-25, self.image_height+25])
            self._styled_polygons.append({"polygon":placeholder_poly, "facecolor":"None", "shade":0.5, "edgecolor":'black', "linewidth":1})
            plot_polygon(placeholder_poly, facecolor="#eeeeee", scale=1.0, axes=axs, edgecolor='black', linewidth=1, zorder_mod=-1)
            self.image_height = 25 + placeholder

        # main plotting loop
        for i, p in enumerate(self._polygons):
            if smoothing:
                poly_to_draw = smooth_polygon(p[1], level=3)
            else:
                poly_to_draw = p[1]

            # look up color for polygon
            if self.outline_by == "residue":
                key_for_color = p[0].get(color_residues_by)
            else:
                key_for_color = p[0].get(self.outline_by)
            fc = color_map.get(key_for_color, sequential_colors[0])
            base_fc = fc # store original color as well as shading

            shade_value = None
            if depth_shading:
                #fc = shade_from_color(fc, i/len(self._polygons), range=shading_range)
                shade_value = p[0].get("depth", 0.5)
                fc = shade_from_color(fc, shade_value, range=shading_range)
            if depth_lines:
                shade_value = p[0].get("depth", 0.5)
                lw = scale_line_width(shade_value, 0, 0.5)
            else:
                lw = line_width
            plot_polygon(poly_to_draw, facecolor=fc, axes=axs, edgecolor=edge_color, linewidth=lw)
            self._styled_polygons.append({"polygon":poly_to_draw, "facecolor":fc, "edgecolor":edge_color, "linewidth":lw, "shade":shade_value, "base_fc":base_fc})

        if save is not None:
            file_ext = os.path.splitext(save)[1].lower()
            assert file_ext in ['.png','.pdf','.svg','.ps'], "Image file extension not supported"
            plt.savefig(save, dpi=dpi, transparent=True, pad_inches=0, bbox_inches='tight')

        if do_show:
            plt.show()
        else:
            return axs

    def export(self, fname, axes=None):
        """Export a pickle object containing styled polygons than can be combined using ``cellscape scene``"""
        assert(len(self._styled_polygons) > 0)

        data = {'polygons':self._styled_polygons, 'name':self.name, 'width':self.image_width, 'height':self.image_height, 'start':self.start_coord, 'end':self.end_coord, 'bottom':self.bottom_coord, 'top':self.top_coord}

        with open('{}.pickle'.format(fname),'wb') as f:
            pickle.dump(data, f)

        print("Exported polygon data to {}.pickle".format(fname), file=sys.stderr)

def make_cartoon(args):
    """Build a cartoon in one-go. Called when running ``cellscape cartoon``."""

    # accept list of chains for backwards-compatibility
    # convert to string e.g. ABCD for current interface
    # can be an issue if chains have more than one letter
    if len(args.chain) == 1:
        chain = args.chain[0]
    else:
        chain = ''.join(args.chain)

    molecule = Cartoon(args.pdb, chain=chain, model=args.model, uniprot=args.uniprot, view=False)

    # open first line to identify view file
    if args.view is not None:
        with open(args.view) as view_f:
            first_line = view_f.readline()
        if first_line[:8] == 'set_view':
            molecule.load_pymol_view(args.view)
        elif first_line[:5] == 'Model':
            molecule.load_chimera_view(args.view)
        else:
            molecule.load_view_matrix(args.view)
    else:
        # if no view matrix provided just use default PDB orientation for now
        molecule.view_matrix = np.identity(3)

    molecule.outline(args.outline_by, depth=args.depth, radius=args.radius, only_annotated=args.only_annotated, only_ca=args.only_ca, depth_contour_interval=args.depth_contour_interval)
    if args.outline_by == "residue" and args.color_by != "same":
        color_residues_by = args.color_by
    else:
        color_residues_by = None

    if len(args.colors) > 0:
        colors = args.colors
    else:
        colors = None
    molecule.plot(do_show=False, axes_labels=args.axes, colors=colors, color_residues_by=color_residues_by, dpi=args.dpi, save=args.save, depth_shading=args.depth_shading, depth_lines=args.depth_lines, edge_color=args.edge_color, line_width=args.line_width)

    if args.export:
        molecule.export(os.path.splitext(args.save)[0])
