import numpy as np
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

def transform_from_nglview(m):
    camera_matrix = np.array(m).reshape(4,4)
    return camera_matrix[:3,:3]/np.linalg.norm(camera_matrix[:3,:3], axis=1), camera_matrix[3,:3]

def group_by(obj, attr):
    d = dict()
    for o in obj:
        a = o[attr]
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

def safe_union(polys):
    # union of list of shapely polygons
    # maybe slower than cascaded_union/unary_union but catches topology errors
    u = polys[0]
    for p in polys:
        try:
            u = u.union(p)
        except:
            pass
    return u

def get_sequential_colors(colors='Set1', n=1):
    cmap = cm.get_cmap(colors)
    sequential_colors = [cmap(x) for x in range(n)]
    #sequential_colors = [cmap(x) for x in np.linspace(0.0,1.0, n)]
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

def plot_polygon(poly, fc='orange', ec='k', linewidth=1, scale=1.0):
    axs = plt.gca()
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
        # load structure with pytraj
        #self.traj = pt.load(file)
        #self.xyz = self.traj.xyz[0]
        #self.view = nv.show_pytraj(self.traj)

        # load structure with biopython
        parser = PDBParser()
        self.structure = parser.get_structure(file, file)[model]
        _all_chains = [c.id for c in self.structure.get_chains()]

        if lc(chain) == "all":
            self.chains = _all_chains
        else:
            self.chains = list(chain)
            for c in _all_chains:
                if c not in self.chains:
                    self.structure.detach_child(c)

        self.view_matrix = []
        if view:
            if 'nglview' not in sys.modules:
                import nglview as nv
            self._structure_to_view = self.structure
            self.view = nv.show_biopython(self._structure_to_view, sync_camera=True)
            self.view._set_sync_camera([self.view])
            self._reflect_y = np.array([[-1,0,0],[0,1,0],[0,0,-1]])

        # uniprot information
        self._uniprot_xml = uniprot

        # leaving in for testing
        self.xyz = np.concatenate([np.array([list(atom.get_vector()) for atom in self.structure[chain].get_atoms()]) for chain in self.chains])

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
                'coord':(i-len(xyz),i),
                'xyz':xyz, # leaving in for testing
                'com':np.mean(xyz, axis=0) # leaving in for testing
                }
        self.coord = np.array(self.coord)
        self.residues_flat = [self.residues[c][i] for c in self.residues for i in self.residues[c]]

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
        nglv_matrix = np.identity(4)
        nglv_matrix[:3,:3] = np.dot(self.view_matrix, self._reflect_y)
        self.view._set_camera_orientation(list(nglv_matrix.flatten()))
        self.view.center()

        #print("view_matrix:\n", self.view_matrix)
        #print("nglv_matrix:\n", nglv_matrix)
        #print("set_camera_orientation:\n", list(nglv_matrix.flatten()))
        #print("_camera_orientation:\n", np.array(self.view._camera_orientation).reshape(4,4))

    def save_view_matrix(self, p):
        self._update_view_matrix()
        np.savetxt(p, self.view_matrix)

    def load_view_matrix(self, p):
        self.view_matrix = np.loadtxt(p)

    def outline(self, by="all", color=None, occlude=True, only_annotated=False):
        self._update_view_matrix()
        self._polygons = []

        self._polygon = so.unary_union([sg.Point(i).buffer(1.5) for i in self._rotate_and_project(self.view_matrix)])

        if by == 'all':
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
                            space_filling = space_filling.union(sg.Point(i).buffer(1.5))
                        #space_filling = so.cascaded_union([sg.Point(i*10+np.random.random(3)).buffer(15) for i in coords])
                        if not view_object:
                            view_object = space_filling
                        else:
                            if view_object.disjoint(space_filling):
                                view_object = view_object.union(space_filling)
                                region_polygons[res.get(by)].append(space_filling)
                            elif view_object.contains(space_filling):
                                pass
                            else:
                                region_polygons[res.get(by)].append(space_filling.difference(view_object))
                                view_object = view_object.union(space_filling)

                for i, (region_name, region_polygon) in enumerate(region_polygons.items()):
                    if not only_annotated or region_name is not None:
                        self._polygons.append(safe_union(region_polygon))
            else:
                residue_groups = group_by(self.residues_flat, by)
                for group_i, (group_name, group_res) in enumerate(sorted(residue_groups.items(), key=not_none)):
                    if not only_annotated or group_name is not None:
                        group_coords = np.concatenate([rotated_coord[range(*res['coord'])] for r in group_res])
                        self._polygons.append(space_filling = so.cascaded_union([sg.Point(i).buffer(1.5) for i in group_coords]))

        print("Outlined some atoms!", file=sys.stdout)

    def plot(self, axes=True, colors=None, smoothing=False):

        # fire up a pyplot
        fig, axs = plt.subplots()
        axs.set_aspect('equal')

        if axes:
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

        # this should all be cleaned up, add more checking
        num_colors_needed = len(self._polygons)
        default_color = 'dodgerblue'
        default_cmap = 'Set1'

        if colors is None:
            if num_colors_needed == 1:
                sequential_colors = [default_color]
            else:
                sequential_colors = get_sequential_colors(colors=default_cmap, n=num_colors_needed)
        else:
            if num_colors_needed == 1:
                sequential_colors = [colors]
            else:
                sequential_colors = get_sequential_colors(colors=colors, n=num_colors_needed)

        for i, p in enumerate(self._polygons):
            if smoothing:
                poly_to_draw = smooth_polygon(p, level=1)
            else:
                poly_to_draw = p
            plot_polygon(poly_to_draw, fc=sequential_colors[i], scale=1.0)

        plt.show()
