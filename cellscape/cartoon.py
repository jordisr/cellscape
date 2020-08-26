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
import re
import os
import sys
import operator
import colorsys
import warnings
from Bio.PDB import *
from Bio.PDB.Polypeptide import three_to_one
from scipy import signal, interpolate
from scipy.spatial.distance import pdist, squareform
import time

from .parse_uniprot_xml import parse_xml, download_uniprot_record
from .parse_alignment import align_pair, overlap_from_alignment
from .util import *

# silence warnings from Biopython that might pop up when loading the PDB
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

def matrix_from_nglview(m):
    """Take flattened 4x4 view matrix from NGLView and convert to 3x3 rotation matrix."""
    camera_matrix = np.array(m).reshape(4,4)
    return camera_matrix[:3,:3]/np.linalg.norm(camera_matrix[:3,:3], axis=1), camera_matrix[3,:3]

def matrix_to_nglview(m):
    """Take 3x3 rotation matrix and convert to flattened 4x4 view matrix for NGLView."""
    nglv_matrix = np.identity(4)
    nglv_matrix[:3,:3] = np.dot(m, np.array([[-1,0,0],[0,1,0],[0,0,-1]]))
    return list(nglv_matrix.flatten())

def group_by(l, key):
    """Take a list of dictionaries and group them according to a key."""
    d = dict()
    for i in l:
        k = key(i)
        if k in d:
            d[k].append(i)
        else:
            d[k] = [i]
    return d

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
    else:
        return p

def ring_coding(ob):
    # https://sgillies.net/2010/04/06/painting-punctured-polygons-with-matplotlib.html
    # The codes will be all "LINETO" commands, except for "MOVETO"s at the
    # beginning of each subpath
    n = len(ob.coords)
    codes = np.ones(n, dtype=Path.code_type) * Path.LINETO
    codes[0] = Path.MOVETO
    return codes

def transform_coord(xy, offset=np.array([0,0]), scale=1.0, flip=False, recenter=None):
    # 2d coordinates
    xy_ = xy
    if recenter is not None:
        # optionally shift coordinates before rotation
        xy_ -= recenter
    if flip:
        xy_ = np.dot(xy_, np.array([[-1,0],[0,-1]]))
        offset_x = np.min(xy_[:,0])
        offset_y = np.min(xy_[:,1])
        xy_ -= np.array([offset_x, offset_y])
    return (xy_+offset)*scale

def polygon_to_path(polygon, min_interior_length=40, offset=np.array([0,0]), scale=1.0, flip=False, recenter=None):
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
    transformed_vertices = transform_coord(vertices, offset=offset, scale=scale, flip=flip, recenter=recenter)
    return Path(transformed_vertices, codes)

def plot_polygon(poly, facecolor='orange', edgecolor='k', linewidth=0.7, axes=None, zorder_mod=0, offset=np.array([0,0]), scale=1.0, flip=False, recenter=None, min_area=7):
    """Draw a Shapely polygon using matplotlib Patches."""
    if axes is None:
        axs = plt.gca()
        axs.set_aspect('equal')
    else:
        axs = axes
    if isinstance(poly, sg.polygon.Polygon):
        if poly.area > min_area:
            path = polygon_to_path(poly, offset=offset, scale=scale, flip=flip, recenter=recenter)
            patch = PathPatch(path, facecolor=facecolor, edgecolor='black', linewidth=linewidth, zorder=3+zorder_mod)
            axs.add_patch(patch)
    elif isinstance(poly, sg.multipolygon.MultiPolygon):
        for p in poly:
            plot_polygon(p, axes=axs, facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, scale=scale, zorder_mod=zorder_mod, offset=offset, flip=flip, recenter=recenter)

def orientation_from_topology(topologies):
    """Infer protein vertical orientation (N->C or C->N) from UniProt topology annotation."""
    first_ex_flag = True
    first_ex = None
    first_cy_flag = True
    first_cy = None
    first_he_flag = True
    first_he = None

    for row in topologies:
        (description, start, end) = row

        if description == 'Extracellular' and first_ex_flag:
            first_ex = (start, end)
            first_ex_flag = False
        elif description == 'Helical' and first_he_flag:
            first_he = (start, end)
            first_he_flag = False
        elif description == 'Cytoplasmic' and first_cy_flag:
            first_cy = (start, end)
            first_cy_flag = False

    # rough heuristic for now, works for single pass transmembrane proteins
    orient_from_topo = 1
    if first_ex is not None and first_cy is not None:
        if first_ex[0] > first_cy[0]:
            orient_from_topo = 1
        elif first_ex[0] < first_cy[0]:
            orient_from_topo = -1

    print("guessed orientation:{}".format(orient_from_topo))
    return(orient_from_topo)

def depth_slices_from_coord(xyz, width):
    # split single xyz Nx3 matrix into list of Nx3 matrices
    binned = (xyz[:,-1]/width).astype(int)
    binned_shifted = binned - np.min(binned)
    num_bins = np.max(binned_shifted)+1

    total_coords = 0
    slice_coords = []

    for i in range(num_bins):
        bin_coords = xyz[binned_shifted == i]
        slice_coords.append(bin_coords)
        total_coords += len(bin_coords)

    assert(len(xyz) == total_coords)
    return slice_coords

def split_on_labels(m, labels):
    num_bins = np.max(labels)+1
    total_coords = 0
    coords = []
    for i in range(num_bins):
        group_coords = m[labels == i]
        coords.append(group_coords)
        total_coords += len(group_coords)
    assert(len(m) == total_coords)
    return coords

def get_z_slice_labels(xyz, width):
    # Take an Nx3 coordinate matrix and return Z bin
    binned = (xyz[:,-1]/width).astype(int)
    return binned - np.min(binned)

class Cartoon:
    """Main object used to build a molecular cartoon from a PDB structure."""
    def __init__(self, file, name=None, model=0, chain="all", uniprot=None, view=True, is_opm=False):

        # descriptive name for the protein
        self.name = name

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

        # assumes PDB is oriented as described here:
        # https://opm.phar.umich.edu/about#features
        self.is_opm = is_opm

        # view matrix and NGLView options
        self.use_nglview = view
        self.view_matrix = []
        if self.use_nglview:
            if 'nglview' not in sys.modules or 'nv' not in sys.modules:
                import nglview as nv
            self._structure_to_view = self.structure
            initial_repr = [
                {"type": "spacefill", "params": {
                    "sele": "protein", "color": "element"
                }}
            ]
            self.view = nv.show_biopython(self._structure_to_view, sync_camera=True, representations=initial_repr)
            self.view.camera = 'orthographic'
            self.view._set_sync_camera([self.view])
            self._reflect_y = np.array([[-1,0,0],[0,1,0],[0,0,-1]])

        # data structure holding residue information
        self.residues = dict()
        self.sequence = dict()
        self.coord = []
        self.ca_atoms = []
        all_atoms = 0
        for chain in self.chains:
            self.sequence[chain] = ""
            self.residues[chain] = dict()
            for res in self.structure[chain]:
                res_id = res.get_full_id()[3][1]
                if res.get_full_id()[3][0][0] == "H": # skip hetatm records
                    continue
                if res.get_resname() not in amino_acid_3letter:
                    continue
                res_aa = amino_acid_3letter[res.get_resname()]
                self.sequence[chain] += res_aa
                residue_atoms = 0
                these_atoms = []
                for a in res:
                    self.coord.append(list(a.get_vector()))
                    these_atoms.append(a.id) # tracking atom identities for now
                    if a.id == "CA":
                        this_ca_atom = all_atoms
                        self.ca_atoms.append(this_ca_atom)
                    all_atoms += 1
                    residue_atoms += 1
                self.residues[chain][res_id] = {
                'chain':chain,
                'id':res_id,
                'amino_acid':res_aa,
                'object':res,
                'coord':(all_atoms-residue_atoms, all_atoms),
                'coord_ca':(this_ca_atom, this_ca_atom+1),
                'atoms':np.array(these_atoms)
                }
        self.coord = np.array(self.coord)
        self.ca_atoms = np.array(self.ca_atoms).astype(int)

        # uniprot information
        if uniprot is not None:
            if os.path.exists(uniprot):
                self._uniprot_xml = uniprot
            elif re.fullmatch(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', uniprot):
                # if file doesn't exist, check it is a valid UniProt ID and download from server
                # using regex from https://www.uniprot.org/help/accession_numbers
                try:
                    self._uniprot_xml = download_uniprot_record(uniprot, "xml", os.getcwd())
                except:
                    sys.exit("Couldn't download UniProt file")
            else:
                self._uniprot_xml = None
                sys.exit("Must specify either a UniProt XML file or a valid UniProt ID")
        else:
            self._uniprot_xml = None

        if self._uniprot_xml is not None:
            self._preprocess_uniprot(self._uniprot_xml)

    def _preprocess_uniprot(self, xml):
        # TODO support more than one XML file (e.g. for different chains?)
        self._uniprot = parse_xml(xml)[0]

        # align PDB and UniProt sequences to find offset
        uniprot_chain = self.chains[0]
        pdb_seq = self.sequence[uniprot_chain]
        uniprot_seq = self._uniprot.sequence
        first_residue_id = sorted(self.residues[uniprot_chain])[0]
        # alignment coordinates are 0-indexed (but PDB numbering and Uniprot ranges are 1-indexed)
        self._uniprot_overlap = np.array(overlap_from_alignment(align_pair(uniprot_seq, pdb_seq))) + np.array([1,1,1,1])
        self._uniprot_offset = self._uniprot_overlap[0] - first_residue_id

        if len(self._uniprot.domains) > 0:
            self._annotate_residues_from_uniprot(self._uniprot.domains, name_key="domain", residues=self.residues[uniprot_chain], offset=self._uniprot_offset)

        if len(self._uniprot.topology) > 0:
            self._annotate_residues_from_uniprot(self._uniprot.topology, name_key="topology", residues=self.residues[uniprot_chain], offset=self._uniprot_offset)

    def _annotate_residues_from_uniprot(self, ranges, name_key, residues, offset=0):
        # pdb_number - offset = up_number
        for row in ranges:
            (name, start, end) = row
            for r in range(start, end+1):
                if (r-offset) in residues:
                    residues[r-offset][name_key] = name

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

    def auto_view(self, n_atoms=100, c_atoms=100, flip=None):
        if flip is None:
            if self._uniprot_xml and hasattr(self._uniprot, "topology"):
                nc_orient = orientation_from_topology(self._uniprot.topology)
            else:
                nc_orient = True
        elif isinstance(flip, bool):
            nc_orient = flip

        # rotate structure so N-C vector is aligned with the vertical axis
        com = np.mean(self.coord, axis=0)
        atoms_ = self.coord - com
        v1 = np.mean(atoms_[:n_atoms], axis=0) - np.mean(atoms_[-c_atoms:], axis=0)
        if nc_orient:
            first_rotation = rotmat(vectors.Vector(v1), vectors.Vector(np.array([0,1,0]))).T
        else:
            first_rotation = rotmat(vectors.Vector(v1), vectors.Vector(np.array([0,-1,0]))).T

        # rotate around Y axis so X axis aligns with longest distance in XZ plane
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
            # waiting a little bit seems to fix it (maybe an issue with sync/refresh rate)
            #self.view.control.orient(nglv_matrix)
            #self.view._camera_orientation = nglv_matrix
            time.sleep(0.5)
            self.view.center()
            #print("After", self.view._camera_orientation)

    def _apply_view_matrix(self):
        # transform atomic coordinates using view matrix
        self.rotated_coord = np.dot(self.coord, self.view_matrix)

    def load_pymol_view(self, file):
        """Read rotation matrix from output of PyMol ``get_view`` command."""
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
        """Read rotation matrix from output of Chimera ``matrixset`` command."""
        matrix = []
        with open(file,'r') as view:
            for line in view.readlines()[1:4]:
                matrix.append(line.split())

        # transpose and remove translation vector
        self.view_matrix = np.array(matrix).astype(float).T[:3]

        # rotate camera in nglview
        self._set_nglview_orientation(self.view_matrix)

    def save_view_matrix(self, p):
        """Save rotation matrix to a NumPy text file."""
        self._update_view_matrix()
        np.savetxt(p, self.view_matrix)

    def load_view_matrix(self, p):
        """Load rotation matrix from a NumPy text file."""
        self.view_matrix = np.loadtxt(p)
        self._set_nglview_orientation(self.view_matrix)

    def set_view_matrix(self, m):
        """Manually set view matrix (3x3)."""
        # TODO abstract other loading functions to use this (i.e. only call _set_nglview_orientation in this function)
        assert m.shape == (3,3)
        self.view_matrix = m
        self._set_nglview_orientation(self.view_matrix)

    def outline(self, by="all", depth=None, depth_contour_interval=3, only_ca=False, only_annotated=False, radius=None, back_outline=False, align_transmembrane=False):
        """Create 2D projection from coordinates and outline atoms."""

        # check options
        assert by in ["all", "residue", "chain", "domain", "topology"], "Option not recognized"
        assert depth in [None, "flat", "contours"], "Option not recognized"
        # depth option doesn't affect by="residues"

        # collapse chain hierarchy into flat list
        self.residues_flat = [self.residues[c][i] for c in self.residues for i in self.residues[c]]

        if self.is_opm:
            self.set_view_matrix(np.array([[1,0,0],[0,0,1],[0,1,0]]))
        elif self.use_nglview:
            self._update_view_matrix()

        # transform atomic coordinates using view matrix
        self._apply_view_matrix()

        # recenter coordinates on lower left edge of bounding box
        offset_x = np.min(self.rotated_coord[:,0])
        if self.is_opm:
            offset_y = 0 # since OPM already aligned to membrane
        else:
            offset_y = np.min(self.rotated_coord[:,1])
        self.rotated_coord -= np.array([offset_x, offset_y, 0])

        # calculate vertical offset for transmembrane proteins
        if self._uniprot_xml and align_transmembrane:
            tm_coordinates = []
            for res in self.residues_flat:
                if res.get("topology","") == "Helical":
                    tm_coordinates.append(np.array(self.rotated_coord[range(*res['coord_ca'])]))
            if len(tm_coordinates) > 0:
                tm_coordinates = np.concatenate(np.array(tm_coordinates))
                tm_com_y = np.mean(tm_coordinates[:,1])
                print("shifted for transmembrane region by {} angstroms".format(tm_com_y))
                self.rotated_coord -= np.array([0, tm_com_y, 0])

        self._rescale_z = lambda z: (z-np.min(self.rotated_coord[:,-1]))/(np.max(self.rotated_coord[:,-1])-np.min(self.rotated_coord[:,-1]))
        self._polygons = []
        self._styled_polygons = []

        # default radius for rendering atoms
        if only_ca and radius is None:
            radius_ = 5
        elif radius is None:
            radius_ = 1.5
        else:
            radius_ = radius

        if by == 'all':
            # space-filling outline of entire molecule
            self.num_groups = 1
            if only_ca:
                coord_to_outline = self.rotated_coord[self.ca_atoms]
            else:
                coord_to_outline = self.rotated_coord
            if depth == "contours":
                slice_coords = split_on_labels(coord_to_outline, get_z_slice_labels(coord_to_outline, width=depth_contour_interval))
                for slice in slice_coords:
                    slice_depth = self._rescale_z(np.mean(slice[:,-1]))
                    self._polygons.append(({"depth":slice_depth}, so.unary_union([sg.Point(i).buffer(radius_) for i in slice])))
            else:
                # depth=None and depth=flat are equivalent for by="all"
                self._polygons.append(({}, so.unary_union([sg.Point(i).buffer(radius_) for i in coord_to_outline])))
        else:
            for res in self.residues_flat:
                # pick range of atomic coordinates out of main data structure
                if only_ca:
                    res_coords = np.array(self.rotated_coord[range(*res['coord_ca'])])
                else:
                    res_coords = np.array(self.rotated_coord[range(*res['coord'])])
                res["xyz"] = res_coords

        if by == 'residue':
            for res in sorted(self.residues_flat, key=lambda res: np.mean(res["xyz"][:,-1])):
                group_outline = so.cascaded_union([sg.Point(i).buffer(radius_) for i in res["xyz"] ])
                res["polygon"] = group_outline
                res["depth"] = self._rescale_z(np.mean(res["xyz"][:,-1]))
                self._polygons.append((res, group_outline))

        elif by in ['domain', 'topology', 'chain']:

            if by in ['domain', 'topology']:
                assert(self._uniprot_xml is not None)

            # TODO comment code and be consistent with variable names group vs region
            residue_groups = group_by(self.residues_flat, key=lambda x: x.get(by))
            self.groups = sorted(residue_groups.keys(), key=lambda x: (x is None, x))

            self.num_groups = len(residue_groups)
            region_atoms = dict() # residue group to atomic indices
            total_atoms = 0
            for k,v in residue_groups.items():
                region_atoms[k] = []
                for res in v:
                    if only_ca:
                        region_atoms[k].extend(range(*res['coord_ca']))
                    else:
                        region_atoms[k].extend(range(*res['coord']))
                region_atoms[k] = np.array(region_atoms[k], dtype=int)
                total_atoms += len(region_atoms[k])
            #print("ATOMS", total_atoms, len(self.rotated_coord)) # for debugging

            if depth is not None:

                slice_labels = get_z_slice_labels(self.rotated_coord, width=depth_contour_interval)
                num_slices = np.max(slice_labels)+1

                if depth == "contours":
                    for s in range(num_slices):
                        for group_i, (group_name, group_res) in enumerate(sorted(residue_groups.items(), key=lambda x: (x[0] is None, x))):
                            if not only_annotated or group_name is not None:
                                atom_indices = region_atoms[group_name]
                                slice_coords = self.rotated_coord[atom_indices][slice_labels[atom_indices] == s]
                                if len(slice_coords) > 0:
                                    slice_depth = self._rescale_z(np.mean(slice_coords[:,-1]))
                                    self._polygons.append(({by:group_name, "depth":slice_depth}, so.unary_union([sg.Point(c).buffer(radius_) for c in slice_coords])))

                elif depth == "flat":
                    empty_polygon = sg.Point((0,0)).buffer(0)
                    view_object = empty_polygon
                    region_polygons = dict()
                    for slice in range(num_slices, 0, -1):
                        for group_i, (group_name, group_res) in enumerate(sorted(residue_groups.items(), key=lambda x: (x[0] is None, x))):
                            if not only_annotated or group_name is not None:
                                atom_indices = region_atoms[group_name]
                                slice_coords = self.rotated_coord[atom_indices][slice_labels[atom_indices] == slice]
                                poly = so.unary_union([sg.Point(c).buffer(radius_) for c in slice_coords])
                                this_difference = poly.difference(view_object)
                                region_polygons[group_name] = region_polygons.get(group_name, empty_polygon).union(this_difference.buffer(0.01))
                                view_object = view_object.union(this_difference.buffer(0.01))

                    for k,v in region_polygons.items():
                        self._polygons.append(({by:k}, v))

            else:
                for group_i, (group_name, group_res) in enumerate(residue_groups.items()):
                    if not only_annotated or group_name is not None:
                        group_coords = self.rotated_coord[region_atoms[group_name]]
                        self._polygons.append(({by:group_name}, so.unary_union([sg.Point(i).buffer(radius_) for i in group_coords])))

        if back_outline:
            self._back_outline =  so.unary_union([p[1].buffer(0.01) for p in self._polygons])
        else:
            self._back_outline = None

        self.outline_by = by
        print("Outlined {} polygons!".format(len(self._polygons)), file=sys.stderr)

    def plot(self, colors=None, axes_labels=False, color_residues_by=None, edge_color="black", line_width=0.7,
        depth_shading=False, depth_lines=False, shading_range=0.4, smoothing=False, do_show=True, axes=None, save=None, dpi=300):
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
            axs.set_aspect('equal')

            if axes_labels:
                axs.xaxis.grid(False)
                axs.yaxis.grid(True)
                axs.axes.xaxis.set_ticklabels([])
                plt.axis('on')
                plt.margins(0,0) # needed to scale axes appropriately?
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
            if self.num_groups == 1:
                sequential_colors = [default_color]
            else:
                sequential_colors = get_sequential_colors(colors=default_cmap, n=self.num_groups)
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
            plot_polygon(self._back_outline, facecolor="None", scale=1.0, axes=axs, edgecolor=edge_color, linewidth=2, zorder_mod=-1)
            self._styled_polygons.append({"polygon":self._back_outline, "facecolor":"None", "edgecolor":edge_color, "linewidth":1})

        # main plotting loop
        for i, p in enumerate(self._polygons):
            if smoothing:
                poly_to_draw = smooth_polygon(p[1], level=1)
            else:
                poly_to_draw = p[1]

            # look up color for polygon
            if self.outline_by == "residue":
                key_for_color = p[0].get(color_residues_by)
            else:
                key_for_color = p[0].get(self.outline_by)
            fc = color_map.get(key_for_color, sequential_colors[0])

            if depth_shading:
                #fc = shade_from_color(fc, i/len(self._polygons), range=shading_range)
                fc = shade_from_color(fc, p[0].get("depth", 0.5), range=shading_range)
            if depth_lines:
                lw = scale_line_width(p[0].get("depth", 0.5), 0, 0.5)
            else:
                lw = line_width
            plot_polygon(poly_to_draw, facecolor=fc, axes=axs, edgecolor=edge_color, linewidth=lw)
            # TODO instead of separate variable, just add style info to polygon?
            self._styled_polygons.append({"polygon":poly_to_draw, "facecolor":fc, "edgecolor":edge_color, "linewidth":lw})

        if save is not None:
            file_ext = os.path.splitext(save)[1].lower()
            assert file_ext in ['.png','.pdf','.svg','.ps'], "Image file extension not supported"
            plt.savefig(save, dpi=dpi, transparent=True, pad_inches=0, bbox_inches='tight')

        if do_show:
            plt.show()
        else:
            return axs

    def export(self, fname, axes=None):
        """Export a pickle object containing styled polygons than can be combined using ``scene``"""
        assert(len(self._styled_polygons) > 0)

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

        data = {'polygons':self._styled_polygons, 'name':self.name, 'width':image_width, 'height':image_height, 'start':start_coord, 'end':end_coord, 'bottom':bottom_coord, 'top':top_coord}

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
