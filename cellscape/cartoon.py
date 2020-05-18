import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib import lines, text, cm
import shapely.geometry as sg
import shapely.ops as so
import pickle
import os
import sys
from Bio.PDB import *
from scipy import signal, interpolate

from .parse_uniprot_xml import parse_xml

def transform_from_nglview(m):
    # take flattened 4x4 matrix from NGLView and convert to 3x3 rotation matrix
    camera_matrix = np.array(m).reshape(4,4)
    return camera_matrix[:3,:3]/np.linalg.norm(camera_matrix[:3,:3], axis=1), camera_matrix[3,:3]

def group_by(obj, attr):
    d = dict()
    for o in obj:
        a = o.get(attr, None)
        if a in d:
            d[a].append(o)
        else:
            d[a] = [o]
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

def get_sequential_colors(colors='Set1', n=1):
    # see uses matplotlib.colors.ColorMap
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

        # load structure with biopython
        parser = PDBParser()
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
            if 'nglview' not in sys.modules:
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

    def _rotate_and_project(self, m):
        return np.dot(self.coord, m)[:,:2]

    def _coord(self, t):
        return self.coord[t[0]:t[1]]

    def _update_view_matrix(self):
        # check if camera orientation has been specified from nglview
        if len(self.view._camera_orientation) == 16:
            m, t = transform_from_nglview(self.view._camera_orientation)
            self.view_matrix = np.dot(m, self._reflect_y)
        elif len(self.view_matrix) == 0:
            self.view_matrix = np.identity(3)

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
        if self.use_nglview:
            nglv_matrix = np.identity(4)
            nglv_matrix[:3,:3] = np.dot(self.view_matrix, self._reflect_y)
            self.view._set_camera_orientation(list(nglv_matrix.flatten()))
            self.view.center()

    def save_view_matrix(self, p):
        self._update_view_matrix()
        np.savetxt(p, self.view_matrix)

    def load_view_matrix(self, p):
        self.view_matrix = np.loadtxt(p)
        # TODO do you need to update NGL camera view here?

    def outline(self, by="all", color=None, occlude=True, only_annotated=False):

        # collapse chain hierarchy into flat list
        self.residues_flat = [self.residues[c][i] for c in self.residues for i in self.residues[c]]

        if self.use_nglview:
            self._update_view_matrix()

        self._polygons = []

        if by == 'all':
            # space-filling outline of entire molecule
            self._polygon = so.unary_union([sg.Point(i).buffer(1.5) for i in self._rotate_and_project(self.view_matrix)])
            self._polygons.append(self._polygon)

        elif by == 'residue':
            rotated_coord = np.dot(self.coord, self.view_matrix)
            z_sorted_residues = sorted(self.residues_flat, key=lambda res: np.mean(rotated_coord[range(*res['coord'])]))

            for res in z_sorted_residues:
                coords = rotated_coord[range(*res['coord'])]
                group_outline = so.cascaded_union([sg.Point(i).buffer(1.5) for i in coords])
                self._polygons.append(group_outline)

        elif by in ['domain', 'topology', 'chain']:

            if by in ['domain', 'topology']:
                assert(self._uniprot_xml is not None)

            rotated_coord = np.dot(self.coord, self.view_matrix)
            residue_groups = group_by(self.residues_flat, by)
            if occlude:
                region_polygons = {c:[] for c in sorted(residue_groups.keys(), key=not_none)}
                view_object = None
                for res in reversed(self.residues_flat):
                    coords = rotated_coord[range(*res['coord'])]
                    if not only_annotated or res.get(by) is not None:
                        space_filling = sg.Point(coords[0]).buffer(1.5)
                        for i in coords:
                            #space_filling = space_filling.union(sg.Point(i).buffer(1.5))
                            space_filling = safe_union(space_filling, sg.Point(i).buffer(1.5))
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
                        self._polygons.append(safe_union_accumulate(region_polygon))
            else:
                residue_groups = group_by(self.residues_flat, by)
                for group_i, (group_name, group_res) in enumerate(sorted(residue_groups.items(), key=not_none)):
                    if not only_annotated or group_name is not None:
                        group_coords = np.concatenate([rotated_coord[range(*res['coord'])] for r in group_res])
                        self._polygons.append(space_filling = safe_union_accumulate([sg.Point(i).buffer(1.5) for i in group_coords]))

        print("Outlined some atoms!", file=sys.stdout)

    def plot(self, axes_labels=False, colors=None, smoothing=False, do_show=True, axes=None):
        """
        mirroring biopython's phylogeny drawing options
        https://biopython.org/DIST/docs/api/Bio.Phylo._utils-module.html
        can optionally pass a matplotlib Axes instance instead of creating a new one
        if do_show is false then return axes object
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

        # this should all be cleaned up, add more checking
        num_colors_needed = len(self._polygons)
        # TODO What if you're repeating colors (e.g. same domains)
        default_color = 'dodgerblue'
        default_cmap = 'Set1'

        if colors is None:
            if num_colors_needed == 1:
                sequential_colors = [default_color]
            else:
                sequential_colors = get_sequential_colors(colors=default_cmap, n=num_colors_needed)
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
                else:
                    if len(colors) > num_colors_needed:
                        # TODO if too many colors, just take the needed ones
                        pass
                    elif len(colors) < num_colors_needed:
                        # TODO if not enough colors, just repeat
                        pass

            elif isinstance(colors, dict):
                pass

        for i, p in enumerate(self._polygons):
            if smoothing:
                poly_to_draw = smooth_polygon(p, level=1)
            else:
                poly_to_draw = p
            plot_polygon(poly_to_draw, fc=sequential_colors[i], scale=1.0, axes=axs)

        if do_show:
            plt.show()
        else:
            return axs

def make_cartoon(args):
    """
    minimal functionality of pdb2svg.py
    doesn't support all arguments
    """

    # accept list of chains for backwards-compatibility
    # convert to string e.g. ABCD for current interface
    if len(args.chain) == 1:
        chain = args.chain[0]
    else:
        chain = ''.join(args.chain)

    molecule = Cartoon(args.pdb, chain=chain, model=args.model, uniprot=args.uniprot, view=False)
    molecule.load_pymol_view(args.view)
    molecule.outline(args.outline_by)
    molecule.plot(do_show=False, linewidth=args.linewidth)
    plt.savefig(args.save+'.'+args.format, axes_labels=args.axes, dpi=args.dpi, transparent=True, pad_inches=0, bbox_inches='tight')
