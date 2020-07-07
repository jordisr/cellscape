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
import colorsys
import warnings
from Bio.PDB import *
from scipy import signal, interpolate
from scipy.spatial.distance import pdist, squareform
import time

from .parse_uniprot_xml import parse_xml

# silence warnings from Biopython that might pop up when loading the PDB
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

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
    exception_counter = 0
    u = polys[0]
    for p in polys:
        try:
            u = u.union(p)
        except:
            exception_counter += 1
            pass
    #print("UNION ACCUMULATE CAUGHT {} EXCEPTION(S)!!".format(exception_counter))
    return u

def safe_union(a, b):
    # catch topology errors
    try:
        c = a.union(b)
    except:
        return a
    return c

def shade_from_color(color, x, range):
    (r, g, b, a) = mcolors.to_rgba(color)
    h, l, s = colorsys.rgb_to_hls(r,g,b)
    l_dark = max(l-range/2, 0)
    l_light = min(l+range/2, 1)
    l_new = l_dark*(1-x) + l_light*x
    return colorsys.hls_to_rgb(h, l_new, s)

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

def plot_polygon(poly, fc='orange', ec='k', linewidth=0.7, scale=1.0, axes=None, zorder_mod=0):
    if axes is None:
        axs = plt.gca()
    else:
        axs = axes
    axs.set_aspect('equal')
    if isinstance(poly, sg.polygon.Polygon):
        xs, ys = poly.exterior.xy
        axs.fill(np.array(xs)/scale, np.array(ys)/scale, alpha=1, fc=fc, ec=ec, linewidth=linewidth, zorder=3+zorder_mod)
    elif isinstance(poly, sg.multipolygon.MultiPolygon):
        for p in poly:
            xs, ys = p.exterior.xy
            axs.fill(np.array(xs)/scale, np.array(ys)/scale, alpha=1, fc=fc, ec=ec, linewidth=linewidth, zorder=3+zorder_mod)
    elif isinstance(poly, (tuple, list)):
        xs, ys = poly
        axs.fill(np.array(xs)/scale, np.array(ys)/scale, alpha=1, fc=fc, ec=ec, linewidth=linewidth, zorder=3+zorder_mod)

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
            initial_repr = [
                {"type": "spacefill", "params": {
                    "sele": "protein", "color": "skyblue"
                }}
            ]
            self.view = nv.show_biopython(self._structure_to_view, sync_camera=True, representations=initial_repr)
            self.view.camera = 'orthographic'
            self.view._set_sync_camera([self.view])
            self._reflect_y = np.array([[-1,0,0],[0,1,0],[0,0,-1]])

        # data structure holding residue information
        self.residues = dict()
        self.coord = []
        self.ca_atoms = []
        all_atoms = 0
        for chain in self.chains:
            self.residues[chain] = dict()
            for res in self.structure[chain]:
                res_id = res.get_full_id()[3][1]
                residue_atoms = 0
                for a in res:
                    self.coord.append(list(a.get_vector()))
                    if a.id == "CA":
                        this_ca_atom = all_atoms
                        self.ca_atoms.append(this_ca_atom)
                    all_atoms += 1
                    residue_atoms += 1
                self.residues[chain][res_id] = {
                'chain':chain,
                'id':res_id,
                'object':res,
                'coord':(all_atoms-residue_atoms, all_atoms),
                'coord_ca':(this_ca_atom, this_ca_atom+1)
                }
        self.coord = np.array(self.coord)
        self.ca_atoms = np.array(self.ca_atoms).astype(int)

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

    # TODO clean up this section
    def align_view(self, v1, v2):
        # rotate structure so v1 is aligned with v2
        r = rotmat(vectors.Vector(v1), vectors.Vector(v2))
        self.view_matrix = r.T
        self._set_nglview_orientation(self.view_matrix)

    def align_view_nc(self, n_atoms=10, c_atoms=10, flip=False):
        # rotate structure so N-C vector is aligned with the vertical axis
        com = np.mean(self.coord, axis=0)
        atoms_ = self.coord - com
        v1 = np.mean(atoms_[:n_atoms], axis=0) - np.mean(atoms_[-c_atoms:], axis=0)
        if not flip:
            self.align_view(v1, np.array([0,1,0]))
        else:
            self.align_view(v1, np.array([0,-1,0]))

    def auto_view(self, n_atoms=100, c_atoms=100, flip=False):
        # rotate structure so N-C vector is aligned with the vertical axis
        com = np.mean(self.coord, axis=0)
        atoms_ = self.coord - com
        v1 = np.mean(atoms_[:n_atoms], axis=0) - np.mean(atoms_[-c_atoms:], axis=0)
        if not flip:
            first_rotation = rotmat(vectors.Vector(v1), vectors.Vector(np.array([0,1,0]))).T
        else:
            first_rotation = rotmat(vectors.Vector(v1), vectors.Vector(np.array([0,-1,0]))).T

        # rotate around Y axis so X axis aligns with longest distance
        rot_coord = np.dot(self.coord, first_rotation)
        com = np.mean(rot_coord, axis=0)
        atoms_ = rot_coord - com
        xz = atoms_[self.ca_atoms][:,[0,2]]
        dist = squareform(pdist(xz))
        max_dist = np.unravel_index(np.argmax(dist, axis=None), dist.shape)
        #print(max_dist, np.max(dist), dist[max_dist[0]][max_dist[1]])
        v2 = atoms_[self.ca_atoms[max_dist[0]]]-atoms_[self.ca_atoms[max_dist[1]]]
        v2[1] = 0
        second_rotation = rotmat(vectors.Vector(v2), vectors.Vector(np.array([1,0,0]))).T

        self.view_matrix = np.dot(first_rotation, second_rotation)
        self._set_nglview_orientation(self.view_matrix)

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

    def _rotate_to_view(self):
        # transform atomic coordinates using view matrix
        self.rotated_coord = np.dot(self.coord, self.view_matrix)

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

    def load_chimera_view(self, file):
        # read rotation matrix from Chimera matrixget command
        matrix = []
        with open(file,'r') as view:
            for line in view.readlines()[1:4]:
                matrix.append(line.split())

        # transpose and remove translation vector
        self.view_matrix = np.array(matrix).astype(float).T[:3]

        # rotate camera in nglview
        self._set_nglview_orientation(self.view_matrix)

    def save_view_matrix(self, p):
        self._update_view_matrix()
        np.savetxt(p, self.view_matrix)

    def load_view_matrix(self, p):
        self.view_matrix = np.loadtxt(p)
        self._set_nglview_orientation(self.view_matrix)

    def outline(self, by="all", color=None, occlude=False, only_ca=False, only_annotated=False, radius=None):

        # collapse chain hierarchy into flat list
        self.residues_flat = [self.residues[c][i] for c in self.residues for i in self.residues[c]]

        #print("Before", self.view._camera_orientation)
        #print("Before2", self.view_matrix)
        if self.use_nglview:
            self._update_view_matrix()
        #print("After", self.view._camera_orientation)
        #print("After2", self.view_matrix)

        # transform atomic coordinates using view matrix
        self._rotate_to_view

        # recenter coordinates on lower left edge of bounding box
        offset_x = np.min(self.rotated_coord[:,0])
        offset_y = np.min(self.rotated_coord[:,1])
        self.rotated_coord -= np.array([offset_x, offset_y, 0])

        self._polygons = []

        # default radius for rendering atoms
        if only_ca and radius is None:
            radius_ = 5
        else:
            radius_ = 1.5

        if by == 'all':
            # space-filling outline of entire molecule
            if only_ca:
                self._polygon = so.unary_union([sg.Point(i).buffer(radius_) for i in self.rotated_coord[self.ca_atoms,:2]])
            else:
                self._polygon = so.unary_union([sg.Point(i).buffer(radius_) for i in self.rotated_coord[:,:2]])
            self._polygons.append(({}, self._polygon))
        else:
            for res in self.residues_flat:
                # pick range of atomic coordinates out of main data structure
                if only_ca:
                    res_coords = np.array(self.rotated_coord[range(*res['coord_ca'])])
                else:
                    res_coords = np.array(self.rotated_coord[range(*res['coord'])])
                res["xyz"] = res_coords

        if by == 'residue':
            for res in self.residues_flat:
                group_outline = so.cascaded_union([sg.Point(i).buffer(radius_) for i in res["xyz"] ])
                res["polygon"] = group_outline
                self._polygons.append((res, group_outline))

        elif by in ['domain', 'topology', 'chain']:

            if by in ['domain', 'topology']:
                assert(self._uniprot_xml is not None)

            residue_groups = group_by(self.residues_flat, key=lambda x: x[by])
            if occlude:
                region_polygons = {c:[] for c in sorted(residue_groups.keys(), key=not_none)}
                view_object = None
                for res in sorted(self.residues_flat, key=lambda res: np.mean(res["xyz"][:,-1]), reverse=True):
                    coords = res["xyz"]
                    if not only_annotated or res.get(by) is not None:
                        space_filling = sg.Point(coords[0]).buffer(radius_)
                        for i in coords:
                            #space_filling = space_filling.union(sg.Point(i).buffer(radius_))
                            space_filling = safe_union(space_filling, sg.Point(i).buffer(radius_))
                        #space_filling = so.cascaded_union([sg.Point(i*10+np.random.random(3)).buffer(15) for i in coords])

                        if not view_object:
                            # initialize if unassigned
                            view_object = space_filling
                        else:
                            if view_object.disjoint(space_filling):
                                # residue doesn't overlap with view so add it to the view
                                view_object = safe_union(view_object, space_filling)
                                region_polygons[res.get(by)].append(space_filling)

                            elif view_object.contains(space_filling):
                                # residue completely covered and not visible
                                pass
                            else:
                                # residue partially occluded
                                #   don't want small holes from other entries
                                # BUG source of TopologyExceptions when accumulating outlines
                                # putting in a threshold for taking the difference seems to work ok
                                #   in one case but I haven't widely tested it
                                difference = space_filling.difference(view_object)
                                #if difference.area > 5:
                                region_polygons[res.get(by)].append(difference.buffer(0.1))
                                view_object = safe_union(view_object, space_filling)

                for i, (region_name, region_polygon) in enumerate(region_polygons.items()):
                    if not only_annotated or region_name is not None:
                        self._polygons.append(({by:region_name}, so.unary_union(region_polygon)))
            else:
                residue_groups = group_by(self.residues_flat, key=lambda x: x[by])
                for group_i, (group_name, group_res) in enumerate(sorted(residue_groups.items(), key=not_none)):
                    if not only_annotated or group_name is not None:
                        group_coords = np.concatenate([self.rotated_coord[range(*r['coord'])] for r in group_res])
                        self._polygons.append(({by:group_name}, so.unary_union([sg.Point(i).buffer(radius_) for i in group_coords])))

        self.outline_by = by
        print("Outlined some atoms!", file=sys.stderr)

    def plot(self, axes_labels=False, colors=None, color_residues_by=None, edge_color="black", line_width=0.7, shading=False, shading_range=0.6, smoothing=False, do_show=True, axes=None, save=None, dpi=300):
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

        if shading:
            z_coord = self.rotated_coord[:,2]
            def rescale_coord(z):
                return (z-np.min(z_coord))/(np.max(z_coord)-np.min(z_coord))

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

        default_color = '#D3D3D3'
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

        if self.outline_by == "residue":
            # TODO clean this section up
            if len(sequential_colors) > 0:
                color_map = {k:sequential_colors[i] for i,k in enumerate(residue_color_groups.keys())}
            else:
                color_map = colors

            z_sorted_residues = sorted(self.residues_flat, key=lambda res: np.mean(res["xyz"][:,-1]))
            for i, p in enumerate(z_sorted_residues):
                if smoothing:
                    poly_to_draw = smooth_polygon(p["polygon"], level=1)
                else:
                    poly_to_draw = p["polygon"]
                color_value = p.get(color_residues_by)
                if isinstance(colors, dict):
                    fc = color_map[color_value]
                else:
                    fc = color_map.get(color_value, sequential_colors[0])

                if shading:
                    fc = shade_from_color(fc, rescale_coord(np.mean(p["xyz"][:,2])), range=shading_range)

                plot_polygon(poly_to_draw, fc=fc, scale=1.0, axes=axs, ec=edge_color, linewidth=line_width)

        else:
            if isinstance(colors, dict):
                sequential_colors = [colors[p[0][self.outline_by]] for p in self._polygons]
            for i, p in enumerate(self._polygons):
                if smoothing:
                    poly_to_draw = smooth_polygon(p[1], level=1)
                else:
                    poly_to_draw = p[1]
                plot_polygon(poly_to_draw, fc=sequential_colors[i], scale=1.0, axes=axs, ec=edge_color, linewidth=line_width)

        if save is not None:
            plt.savefig(save, dpi=dpi, transparent=True, pad_inches=0, bbox_inches='tight')

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

    # open first line to identify view file
    with open(args.view) as view_f:
        first_line = view_f.readline()
    if first_line[:8] == 'set_view':
        molecule.load_pymol_view(args.view)
    elif first_line[:5] == 'Model':
        molecule.load_chimera_view(args.view)
    else:
        molecule.load_view_matrix(args.view)

    molecule.outline(args.outline_by, occlude=args.occlude, radius=args.radius)
    if args.outline_by == "residue" and args.color_by != "same":
        color_residues_by = args.color_by
    else:
        color_residues_by = None

    if len(args.colors) > 0:
        colors = args.colors
    else:
        colors = None
    molecule.plot(do_show=False, axes_labels=args.axes, colors=colors, color_residues_by=color_residues_by, dpi=args.dpi, save="{}.{}".format(args.save, args.format), shading=True, edge_color=args.edge_color, line_width=args.line_width)
