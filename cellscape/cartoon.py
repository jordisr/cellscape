import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib import lines, text, cm
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import shapely.geometry as sg
import shapely.ops as so
import pickle
import os
import sys
import operator
from Bio.PDB import *
from scipy import signal, interpolate
import time

from .parse_uniprot_xml import parse_xml

def matrix_from_nglview(m):
    # take flattened 4x4 matrix from NGLView and convert to 3x3 rotation matrix
    camera_matrix = np.array(m).reshape(4,4)
    return camera_matrix[:3,:3]/np.linalg.norm(camera_matrix[:3,:3], axis=1), camera_matrix[3,:3]

def matrix_to_nglview(m):
    # take 3x3 rotation matrix and convert to flattened 4x4 for NGLView
    nglv_matrix = np.identity(4)
    nglv_matrix[:3,:3] = np.dot(m, np.array([[-1,0,0],[0,1,0],[0,0,-1]]))
    return list(nglv_matrix.flatten())

def group_by(l, key):
    d = dict()
    for i in l:
        k = key(i)
        if k in d:
            d[k].append(i)
        else:
            d[k] = [i]
    return d

def not_none(x):
    if x is None:
        return ""
    else:
        if isinstance(x, str):
            return x
        else:
            if x[0] is None:
                return ""
            else:
                return str(x[1])

def safe_union_accumulate(polys):
    # union of list of shapely polygons
    # maybe slower than cascaded_union/unary_union but catches topology errors
    # TODO figure out a better way of handling these exceptions
    u = polys[0]
    for p in polys:
        try:
            u = u.union(p)
        except:
            pass
    return u

def safe_union(a, b):
    # catch topology errors
    try:
        c = a.union(b)
    except:
        return a
    return c

def shades_from_rgb(color, w=0.3):
    # define custom Lighter Color => Color => Darker Color cmap
    if istype(color, str):
        assert(h[0] == "#" and len(h) == 7)
        r = int(h[1:3], 16)
        g = int(h[3:5], 16)
        b = int(h[5:7], 16)
    elif istype(color, (tuple, list)):
        assert(len(color) == 3)
        (r, g, b) = h
    h, l, s = colorsys.rgb_to_hls(r/255,g/255,b/255)

    # lighter and darker versions of color in HLS space
    c3 = (h, min(l+(1-w)*l, 1), s)
    c2 = (h, l, s)
    c1 = (h, max(l-(1-w)*l, 0), s)

    # convert back to RGB and return colormap
    colors = [colorsys.hls_to_rgb(c[0], c[1], c[2]) for c in [c1, c2, c3]]
    return LinearSegmentedColormap.from_list("shade", colors)

def get_sequential_colors(colors='Set1', n=1):
    # sample n colors from a colormap
    # uses matplotlib.colors.ColorMap.N to distinguish continuous/discrete
    cmap = cm.get_cmap(colors)
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
    else:
        return p

def plot_polygon(poly, fc='orange', ec='k', linewidth=0.7, scale=1.0, axes=None):
    if axes is None:
        axs = plt.gca()
    else:
        axs = axes
    axs.set_aspect('equal')
    if isinstance(poly, sg.polygon.Polygon):
        xs, ys = poly.exterior.xy
        axs.fill(np.array(xs)/scale, np.array(ys)/scale, alpha=1, fc=fc, ec=ec, linewidth=linewidth, zorder=3)
    elif isinstance(poly, sg.multipolygon.MultiPolygon):
        for p in poly:
            xs, ys = p.exterior.xy
            axs.fill(np.array(xs)/scale, np.array(ys)/scale, alpha=1, fc=fc, ec=ec, linewidth=linewidth, zorder=3)
    elif isinstance(poly, (tuple, list)):
        xs, ys = poly
        axs.fill(np.array(xs)/scale, np.array(ys)/scale, alpha=1, fc=fc, ec=ec, linewidth=linewidth, zorder=3)

class Cartoon:
    def __init__(self, file, model=0, chain="all", uniprot=None, view=True):

        # check whether outline has been generated yet
        self.outline_by = None

        # load structure with biopython
        if file[-3:] in ["cif", "mcif"]:
            parser = MMCIFParser()
        elif file[-3:] in ["pdb", "ent"]:
            parser = PDBParser()
        else:
            sys.exit("File format not recognized!")
        self.structure = parser.get_structure(file, file)[model]
        _all_chains = [c.id for c in self.structure.get_chains()]

        # eliminate undesired chains from the biopython object
        if chain.lower() == "all":
            self.chains = _all_chains
        else:
            self.chains = list(chain)
            for c in _all_chains:
                if c not in self.chains:
                    self.structure.detach_child(c)

        # view matrix and NGLView options
        self.use_nglview = view
        self.view_matrix = []
        if self.use_nglview:
            if 'nglview' not in sys.modules or 'nv' not in sys.modules:
                import nglview as nv
            self._structure_to_view = self.structure
            self.view = nv.show_biopython(self._structure_to_view, sync_camera=True)
            self.view._set_sync_camera([self.view])
            self._reflect_y = np.array([[-1,0,0],[0,1,0],[0,0,-1]])

        # data structure holding residue information
        self.residues = dict()
        self.coord = []
        i = 0
        for chain in self.chains:
            self.residues[chain] = dict()
            for res in self.structure[chain]:
                res_id = res.get_full_id()[3][1]
                xyz = [list(a.get_vector()) for a in res] # leaving in for testing
                for coord in xyz:
                    self.coord.append(coord)
                    i += 1
                self.residues[chain][res_id] = {
                'chain':chain,
                'id':res_id,
                'object':res,
                'coord':(i-len(xyz),i)
                }
        self.coord = np.array(self.coord)

        # uniprot information
        self._uniprot_xml = uniprot
        if self._uniprot_xml is not None:
            self._preprocess_uniprot(self._uniprot_xml)

    def _preprocess_uniprot(self, xml):
        # TODO support more than one XML file (e.g. for differnet chains)
        self._uniprot = parse_xml(xml[0])[0]

        # TODO add sequence alignment to find this automatically
        uniprot_chain = self.chains[0]
        uniprot_offset = 0

        if len(self._uniprot.domains) > 0:
            self._annotate_residues_from_uniprot(self._uniprot.domains, name_key="domain", residues=self.residues[uniprot_chain], offset=uniprot_offset)

        if len(self._uniprot.topology) > 0:
            self._annotate_residues_from_uniprot(self._uniprot.topology, name_key="topology", residues=self.residues[uniprot_chain], offset=uniprot_offset)

    def _annotate_residues_from_uniprot(self, ranges, name_key, residues, offset=0):
        for row in ranges:
            (name, start, end) = row
            for r in range(start, end+1):
                if (r+offset) in residues:
                    residues[r+offset][name_key] = name

    def _coord(self, t):
        return self.coord[t[0]:t[1]]

    def _update_view_matrix(self):
        # check if camera orientation has been specified from nglview
        if len(self.view._camera_orientation) == 16:
            m, t = matrix_from_nglview(self.view._camera_orientation)
            self.view_matrix = np.dot(m, self._reflect_y)
        elif len(self.view_matrix) == 0:
            self.view_matrix = np.identity(3)

    def _set_nglview_orientation(self, m):
        # m is 3x3 rotation matrix
        if self.use_nglview:
            nglv_matrix = matrix_to_nglview(m)
            #print("Before", self.view._camera_orientation)
            self.view._set_camera_orientation(nglv_matrix)
            # having a bug where setting camera orientation does nothing
            # waiting a little bit seems to fix it (maybe an issue with sync/refresh)
            #self.view.control.orient(nglv_matrix)
            #self.view._camera_orientation = nglv_matrix
            time.sleep(0.1)
            self.view.center()
            #print("After", self.view._camera_orientation)

    def load_pymol_view(self, file):
        # read rotation matrix from PyMol get_view command
        matrix = []
        with open(file,'r') as view:
            for line in view:
                fields = line.split(',')
                if len(fields) == 4:
                    matrix.append(list(map(float,fields[:3])))
        self.view_matrix = np.array(matrix)[:3]

        # rotate camera in nglview
        self._set_nglview_orientation(self.view_matrix)

    def test(self):
        print(self.view._camera_orientation)

    def save_view_matrix(self, p):
        self._update_view_matrix()
        np.savetxt(p, self.view_matrix)

    def load_view_matrix(self, p):
        self.view_matrix = np.loadtxt(p)
        self._set_nglview_orientation(self.view_matrix)

    def outline(self, by="all", color=None, occlude=True, only_annotated=False, radius=1.5):

        # collapse chain hierarchy into flat list
        self.residues_flat = [self.residues[c][i] for c in self.residues for i in self.residues[c]]

        #print("Before", self.view._camera_orientation)
        #print("Before2", self.view_matrix)
        if self.use_nglview:
            self._update_view_matrix()
        #print("After", self.view._camera_orientation)
        #print("After2", self.view_matrix)

        # transform atomic coordinates using view matrix
        self.rotated_coord = np.dot(self.coord, self.view_matrix)

        # recenter coordinates on lower left edge of bounding box
        offset_x = np.min(self.rotated_coord[:,0])
        offset_y = np.min(self.rotated_coord[:,1])
        self.rotated_coord -= np.array([offset_x, offset_y, 0])

        self._polygons = []

        if by == 'all':
            # space-filling outline of entire molecule
            self._polygon = so.unary_union([sg.Point(i).buffer(radius) for i in self.rotated_coord[:,:2]])
            self._polygons.append(({}, self._polygon))

        elif by == 'residue':
            for res in self.residues_flat:
                # pick range of atomic coordinates out of main data structure
                res_coords = np.array(self.rotated_coord[range(*res['coord'])])
                group_outline = so.cascaded_union([sg.Point(i).buffer(radius) for i in res_coords])
                res["xyz"] = res_coords
                res["polygon"] = group_outline
                self._polygons.append((res, group_outline))

        elif by in ['domain', 'topology', 'chain']:

            if by in ['domain', 'topology']:
                assert(self._uniprot_xml is not None)

            residue_groups = group_by(self.residues_flat, key=lambda x: x[by])
            if occlude:
                region_polygons = {c:[] for c in sorted(residue_groups.keys(), key=not_none)}
                view_object = None
                for res in reversed(self.residues_flat):
                    coords = self.rotated_coord[range(*res['coord'])]
                    if not only_annotated or res.get(by) is not None:
                        space_filling = sg.Point(coords[0]).buffer(radius)
                        for i in coords:
                            #space_filling = space_filling.union(sg.Point(i).buffer(radius))
                            space_filling = safe_union(space_filling, sg.Point(i).buffer(radius))
                        #space_filling = so.cascaded_union([sg.Point(i*10+np.random.random(3)).buffer(15) for i in coords])
                        if not view_object:
                            view_object = space_filling
                        else:
                            if view_object.disjoint(space_filling):
                                #view_object = view_object.union(space_filling)
                                view_object = safe_union(view_object, space_filling)
                                region_polygons[res.get(by)].append(space_filling)
                            elif view_object.contains(space_filling):
                                pass
                            else:
                                region_polygons[res.get(by)].append(space_filling.difference(view_object))
                                #view_object = view_object.union(space_filling)
                                view_object = safe_union(view_object, space_filling)

                for i, (region_name, region_polygon) in enumerate(region_polygons.items()):
                    if not only_annotated or region_name is not None:
                        self._polygons.append(({by:region_name}, safe_union_accumulate(region_polygon)))
            else:
                residue_groups = group_by(self.residues_flat, key=lambda x: x[by])
                for group_i, (group_name, group_res) in enumerate(sorted(residue_groups.items(), key=not_none)):
                    if not only_annotated or group_name is not None:
                        group_coords = np.concatenate([self.rotated_coord[range(*r['coord'])] for r in group_res])
                        self._polygons.append(({by:region_name}, safe_union_accumulate([sg.Point(i).buffer(radius) for i in group_coords])))

        self.outline_by = by
        print("Outlined some atoms!", file=sys.stderr)

    def plot(self, axes_labels=False, colors=None, color_residues_by=None, shading=True, smoothing=False, do_show=True, axes=None, save=None, dpi=300, format="pdf"):
        """
        mirroring biopython's phylogeny drawing options
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

        if axes is None:
            # create a new matplotlib figure if none provided
            fig, axs = plt.subplots()
            axs.set_aspect('equal')

            if axes_labels:
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
        else:
            assert(isinstance(axes, matplotlib.axes.Axes))
            axs = axes
        self._axes= axs

        # TODO this should all be cleaned up, add more checking
        # TODO What if you're repeating colors (e.g. same domains)
        if self.outline_by == "residue":
            if color_residues_by is None:
                num_colors_needed = 1
                residue_color_groups = {"all":self.residues_flat}
            else:
                residue_color_groups = group_by(self.residues_flat, lambda x: x.get(color_residues_by))
                num_colors_needed = len(residue_color_groups)
        else:
            num_colors_needed = len(self._polygons)

        default_color = 'lightgray'
        default_cmap = 'Set1'
        named_colors = [*mcolors.BASE_COLORS.keys(), *mcolors.TABLEAU_COLORS.keys(), *mcolors.CSS4_COLORS.keys(), *mcolors.XKCD_COLORS.keys()]

        if colors is None:
            if num_colors_needed == 1:
                sequential_colors = [default_color]
            else:
                sequential_colors = get_sequential_colors(colors=default_cmap, n=num_colors_needed)
        else:
            if isinstance(colors, dict):
                sequential_colors = []
            else:
                if isinstance(colors, str):
                    if num_colors_needed == 1:
                        sequential_colors = [colors]
                    else:
                        sequential_colors = get_sequential_colors(colors=colors, n=num_colors_needed)
                elif isinstance(colors, (list, tuple)):
                    if num_colors_needed == 1:
                        sequential_colors = [colors[0]]
                    elif num_colors_needed == len(colors):
                        sequential_colors = colors

        if len(sequential_colors) > 0:
            color_map = {k:sequential_colors[i] for i,k in enumerate(residue_color_groups.keys())}

        if self.outline_by == "residue":
            # TODO clean this up
            z_sorted_residues = sorted(self.residues_flat, key=lambda res: np.mean(res["xyz"][:,-1]))
            for i, p in enumerate(z_sorted_residues):
                if smoothing:
                    poly_to_draw = smooth_polygon(p["polygon"], level=1)
                else:
                    poly_to_draw = p["polygon"]
                color_value = p.get(color_residues_by)
                if isinstance(colors, dict):
                    fc = colors[color_value]
                else:
                    fc = color_map.get(color_value, default_color)

                plot_polygon(poly_to_draw, fc=fc, scale=1.0, axes=axs)

        else:
            if isinstance(colors, dict):
                sequential_colors = [colors[p[0][self.outline_by]] for p in self._polygons]
            for i, p in enumerate(self._polygons):
                if smoothing:
                    poly_to_draw = smooth_polygon(p[1], level=1)
                else:
                    poly_to_draw = p[1]
                plot_polygon(poly_to_draw, fc=sequential_colors[i], scale=1.0, axes=axs)

        if save is not None:
            plt.savefig("{}.{}".format(save, format), dpi=dpi, transparent=True, pad_inches=0, bbox_inches='tight')

        if do_show:
            plt.show()
        else:
            return axs

    def export(self, fname, axes=None):
        """
        Find polygons in a matplotlib axes (after self.plot()) and export to pickle for scene building
        """

        if axes is None:
            ax = self._axes
        else:
            ax = axes

        # output summary line
        image_width = np.max(self.rotated_coord[:,0]) - np.min(self.rotated_coord[:,0])
        image_height = np.max(self.rotated_coord[:,1]) - np.min(self.rotated_coord[:,1])
        start_coord = np.mean(self.rotated_coord[:50])
        end_coord = np.mean(self.rotated_coord[:-50])
        bottom_coord = min(self.rotated_coord, key=operator.itemgetter(1))
        top_coord = max(self.rotated_coord, key=operator.itemgetter(1))
        #print('\t'.join(map(str, [image_height, image_width, bottom_coord, top_coord])))

        data = {'polygons':[], 'width':image_width, 'height':image_height, 'start':start_coord, 'end':end_coord, 'bottom':bottom_coord, 'top':top_coord}
        # TODO does this account for patches with holes
        for o in ax.findobj(matplotlib.patches.Polygon):
            data['polygons'].append(o)

        with open('{}.pickle'.format(fname),'wb') as f:
            pickle.dump(data, f)

def make_cartoon(args):
    """
    basic functionality of pdb2svg.py (doesn't support all arguments)
    """

    # accept list of chains for backwards-compatibility
    # convert to string e.g. ABCD for current interface
    if len(args.chain) == 1:
        chain = args.chain[0]
    else:
        chain = ''.join(args.chain)

    molecule = Cartoon(args.pdb, chain=chain, model=args.model, uniprot=args.uniprot, view=False)
    if args.view[-3:] == "npy":
        molecule.load_view_matrix(args.view)
    else:
        molecule.load_pymol_view(args.view)
    molecule.outline(args.outline_by, occlude=args.occlude, radius=args.radius)
    if args.outline_by == "residue":
        color_residues_by = args.color_by
    molecule.plot(do_show=False, axes_labels=args.axes, color_residues_by=color_residues_by, dpi=args.dpi, save=args.save, format=args.format)
