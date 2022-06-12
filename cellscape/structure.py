import numpy as np
import shapely.geometry as sg
import shapely.ops as so
import re
import os
import sys
import operator
import warnings
from Bio.PDB import rotmat, vectors, MMCIFParser, PDBParser
from scipy.spatial.distance import pdist, squareform
import time

import cellscape
from cellscape.util import amino_acid_3letter, group_by
from cellscape.parse_uniprot_xml import parse_xml, download_uniprot_record
from cellscape.parse_alignment import align_pair, overlap_from_alignment

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
    nc_orient = True
    if first_ex is not None and first_cy is not None:
        if first_ex[0] < first_cy[0]:
            nc_orient = True # N->C (top to bottom)
        elif first_ex[0] > first_cy[0]:
            nc_orient = False # C->N (top to bottom)

    return(nc_orient)

def orientation_from_ptm(ptm):
    """Assumes signal peptide is on the cytoplasmic/membrane side with the chain extracellular"""

    nc_orient = True
    if ('chain' in ptm) and ('signal peptide' in ptm):
        if ptm['signal peptide'][0] < ptm['chain'][0]:
            nc_orient = True
        else:
            nc_orient = False

    return(nc_orient)

def depth_slices_from_coord(xyz, width):
    """Split single xyz Nx3 matrix into list of Nx3 matrices"""
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

def get_z_slice_labels(xyz, width):
    """Take an Nx3 coordinate matrix and return Z bin"""
    binned = (xyz[:,-1]/width).astype(int)
    return binned - np.min(binned)

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

def get_dimensions(xy, end_window=50):
    dimensions = {}
    dimensions['width'] = np.max(xy[:,0]) - np.min(xy[:,0])
    dimensions['height'] = np.max(xy[:,1]) - np.min(xy[:,1])
    dimensions['start'] = np.mean(xy[:end_window])
    dimensions['end'] = np.mean(xy[:-end_window])
    dimensions['bottom'] = min(xy, key=operator.itemgetter(1))
    dimensions['top'] = max(xy, key=operator.itemgetter(1))
    return dimensions

class Structure:
    """Load PDB/MMCIF structure and handle NGLView instance"""
    def __init__(self, file, name=None, model=0, chain="all", uniprot=None, view=True, is_opm=False, res_start=None, res_end=None):

        # descriptive name for the protein, otherwise use file
        if name is None:
            self.name = os.path.basename(file)
        else:
            self.name = name

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

        # take chain start and end for first chain
        if res_start is not None and res_end is not None:
            assert(res_end > res_start)
            for res in list(self.structure[_all_chains[0]]):
                res_id = res.get_full_id()[3][1]
                if (res_id < res_start) or (res_id > res_end):
                    self.structure[_all_chains[0]].detach_child(res.get_id())

        # BUG with some biopython structures not loading in nglview
        # can be fixed by resetting disordered flags
        # could this cause problems later on?
        for chain in self.structure:
            for residue in chain:
                for atom in residue.get_unpacked_list():
                    atom.disordered_flag = 0

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
        self.backbone_atoms = []
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
                backbone_atoms = []
                for a in res:
                    self.coord.append(list(a.get_vector()))
                    these_atoms.append(a.id) # tracking atom identities for now
                    if a.id == "CA":
                        this_ca_atom = all_atoms
                        self.ca_atoms.append(this_ca_atom)
                    if a.id in ["CA", "N", "C", "O"]:
                        backbone_atoms.append(all_atoms)
                        self.backbone_atoms.append(all_atoms)
                    all_atoms += 1
                    residue_atoms += 1
                self.residues[chain][res_id] = {
                'chain':chain,
                'id':res_id,
                'amino_acid':res_aa,
                'object':res,
                'coord':(all_atoms-residue_atoms, all_atoms),
                'coord_ca':(this_ca_atom, this_ca_atom+1),
                'coord_backbone':np.array(backbone_atoms),
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

    def _update_view_matrix(self):
        # check if camera orientation has been specified from nglview
        if len(self.view._camera_orientation) == 16:
            m, t = matrix_from_nglview(self.view._camera_orientation)
            self.view_matrix = np.dot(m, self._reflect_y)
        elif len(self.view_matrix) == 0:
            self.view_matrix = np.identity(3)

    def align_view(self, v1, v2):
        # rotate structure so v1 is aligned with v2
        r = rotmat(vectors.Vector(v1), vectors.Vector(v2))
        view_matrix = r.T
        self.set_view_matrix(view_matrix)

    def align_view_nc(self, n_atoms=10, c_atoms=10, flip=False):
        """Rotate structure so N-C vector is aligned with the vertical axis"""
        com = np.mean(self.coord, axis=0)
        atoms_ = self.coord - com
        v1 = np.mean(atoms_[:n_atoms], axis=0) - np.mean(atoms_[-c_atoms:], axis=0)
        if not flip:
            self.align_view(v1, np.array([0,1,0]))
        else:
            self.align_view(v1, np.array([0,-1,0]))

    def auto_view(self, n_atoms=100, c_atoms=100, flip=None):
        """Infer orientation from UniProt data."""
        # TODO should be same as align_view_nc if no UniProt data?
        # TODO abstract with align_view?
        # TODO abstract rotmat to separate function e.g. get_rotation_matrix()
        if flip is None:
            if self._uniprot_xml and len(self._uniprot.topology) > 0:
                print("orienting based on topology...")
                nc_orient = orientation_from_topology(self._uniprot.topology)
            elif self._uniprot_xml and len(self._uniprot.ptm) > 0:
                print("orienting based on ptm...")
                nc_orient = orientation_from_ptm(self._uniprot.ptm)
            else:
                nc_orient = True
        elif isinstance(flip, bool):
            nc_orient = flip
        print("guessed N>C orientation? {}".format(nc_orient))
        self.nc_orient = nc_orient

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

        view_matrix = np.dot(first_rotation, second_rotation)
        self.set_view_matrix(view_matrix)

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
        view_matrix = np.array(matrix)[:3]
        self.set_view_matrix(view_matrix)

    def load_chimera_view(self, file):
        """Read rotation matrix from output of Chimera ``matrixset`` command."""
        matrix = []
        with open(file,'r') as view:
            for line in view.readlines()[1:4]:
                matrix.append(line.split())

        # transpose and remove translation vector
        view_matrix = np.array(matrix).astype(float).T[:3]
        self.set_view_matrix(view_matrix)

    def save_view_matrix(self, p):
        """Save rotation matrix to a NumPy text file."""
        self._update_view_matrix()
        np.savetxt(p, self.view_matrix)

    def load_view_matrix(self, p):
        """Load rotation matrix from a NumPy text file."""
        view_matrix = np.loadtxt(p)
        self.set_view_matrix(view_matrix)

    def set_view_matrix(self, m):
        """Manually set view matrix (3x3)."""
        assert m.shape == (3,3)
        self.view_matrix = m
        self._set_nglview_orientation(self.view_matrix)

    def outline(self, by="all", depth=None, depth_contour_interval=3, only_backbone=False, only_ca=False, only_annotated=False, radius=None, back_outline=False, align_transmembrane=False):
        """Create 2D projection from coordinates and outline atoms."""
        # TODO: expand docstring

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
        polygons = []
        groups = {}
        self._group_outlines = []

        # default radius for rendering atoms
        if only_ca and radius is None:
            radius_ = 5
        elif only_backbone and radius is None:
            radius_ = 4
        elif radius is None:
            radius_ = 1.5
        else:
            radius_ = radius

        if by == 'all':
            # space-filling outline of entire molecule
            self.num_groups = 1
            if only_ca:
                coord_to_outline = self.rotated_coord[self.ca_atoms]
            elif only_backbone:
                coord_to_outline = self.rotated_coord[self.backbone_atoms]
            else:
                coord_to_outline = self.rotated_coord
            if depth == "contours":
                slice_coords = split_on_labels(coord_to_outline, get_z_slice_labels(coord_to_outline, width=depth_contour_interval))
                for slice in slice_coords:
                    slice_depth = self._rescale_z(np.mean(slice[:,-1]))
                    polygons.append(({"depth":slice_depth}, so.unary_union([sg.Point(i).buffer(radius_) for i in slice])))
            else:
                # depth=None and depth=flat are equivalent for by="all"
                polygons.append(({}, so.unary_union([sg.Point(i).buffer(radius_) for i in coord_to_outline])))
        else:
            for res in self.residues_flat:
                # pick range of atomic coordinates out of main data structure
                if only_ca:
                    res_coords = np.array(self.rotated_coord[range(*res['coord_ca'])])
                elif only_backbone:
                    res_coords = np.array(self.rotated_coord[range(*res['coord_backbone'])])
                else:
                    res_coords = np.array(self.rotated_coord[range(*res['coord'])])
                res["xyz"] = res_coords

        if by == 'residue':
            for res in sorted(self.residues_flat, key=lambda res: np.mean(res["xyz"][:,-1])):
                group_outline = so.cascaded_union([sg.Point(i).buffer(radius_) for i in res["xyz"] ])
                res["polygon"] = group_outline
                res["depth"] = self._rescale_z(np.mean(res["xyz"][:,-1]))
                polygons.append((res, group_outline))
            self.num_groups = 1

        elif by in ['domain', 'topology', 'chain']:

            if by in ['domain', 'topology']:
                assert(self._uniprot_xml is not None)

            # TODO comment code and be consistent with variable names group vs region
            residue_groups = group_by(self.residues_flat, key=lambda x: x.get(by))
            groups = sorted(residue_groups.keys(), key=lambda x: (x is None, x))

            self.num_groups = len(residue_groups)
            region_atoms = dict() # residue group to atomic indices
            total_atoms = 0
            for k,v in residue_groups.items():
                region_atoms[k] = []
                for res in v:
                    if only_ca:
                        region_atoms[k].extend(range(*res['coord_ca']))
                    elif only_backbone:
                        region_atoms[k].extend(range(*res['coord_backbone']))
                    else:
                        region_atoms[k].extend(range(*res['coord']))
                region_atoms[k] = np.array(region_atoms[k], dtype=int)
                total_atoms += len(region_atoms[k])

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
                                    slice_outline = so.unary_union([sg.Point(c).buffer(radius_) for c in slice_coords])
                                    polygons.append(({by:group_name, "depth":slice_depth}, slice_outline))

                    # back outline to highlight each group's contours... just duplicating depth==flat code here
                    if back_outline:
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

                        for v in region_polygons.values():
                            self._group_outlines.append(v)

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
                        polygons.append(({by:k}, v))

            else:
                for group_i, (group_name, group_res) in enumerate(residue_groups.items()):
                    if not only_annotated or group_name is not None:
                        group_coords = self.rotated_coord[region_atoms[group_name]]
                        polygons.append(({by:group_name}, so.unary_union([sg.Point(i).buffer(radius_) for i in group_coords])))

        if back_outline:
            self._back_outline =  so.unary_union([p[1].buffer(0.01) for p in polygons])
        else:
            self._back_outline = None

        print("Outlined {} polygons!".format(len(polygons)), file=sys.stderr)

        return cellscape.Cartoon(self.name, polygons, self.residues_flat, by, self._back_outline, self._group_outlines, self.num_groups, get_dimensions(self.rotated_coord), groups)
