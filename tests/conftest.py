"""Shared fixtures: coordinate files in every supported format, built from the bundled PDB."""
import gzip
import shutil
from pathlib import Path

import pytest
from Bio.PDB import MMCIFIO, PDBParser

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"
IG_PDB = EXAMPLES / "ig" / "1igt.pdb"
IG_VIEW = EXAMPLES / "ig" / "view"

# chain IDs used in the multi-letter mmCIF file (1igt has chains A, B, C, D)
MULTI_LETTER_CHAINS = {"A": "AA", "B": "BB", "C": "C", "D": "D"}


def gzip_copy(src, dst):
    with open(src, "rb") as f_in, gzip.open(dst, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    return dst


def write_cif(structure, path):
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(path))
    return path


@pytest.fixture(scope="session")
def coord_files(tmp_path_factory):
    """The 1igt structure written as PDB and mmCIF, in each extension variant we accept."""
    d = tmp_path_factory.mktemp("coords")
    structure = PDBParser(QUIET=True).get_structure("1igt", str(IG_PDB))
    cif = write_cif(structure, d / "1igt.cif")

    files = {
        "pdb": IG_PDB,
        "pdb.gz": gzip_copy(IG_PDB, d / "1igt.pdb.gz"),
        "cif": cif,
        "CIF": shutil.copy(cif, d / "1IGT.CIF"),
        "mmcif": shutil.copy(cif, d / "1igt.mmcif"),
        "cif.gz": gzip_copy(cif, d / "1igt.cif.gz"),
    }

    # same structure with some two-letter chain IDs, which only mmCIF can hold
    for chain in list(structure[0]):
        chain.id = MULTI_LETTER_CHAINS[chain.id]
    files["multi"] = write_cif(structure, d / "1igt_multi.cif")
    return {k: Path(v) for k, v in files.items()}
