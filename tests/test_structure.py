"""Loading coordinate files and selecting chains with cellscape.Structure."""
import numpy as np
import pytest

import cellscape
from conftest import IG_PDB, MULTI_LETTER_CHAINS

# atom counts per chain in 1igt.pdb, as read by cellscape (amino acids only)
ATOMS_PER_CHAIN = {"A": 2033, "B": 4232, "C": 2033, "D": 4232}


def load(path, **kwargs):
    return cellscape.Structure(str(path), view=False, **kwargs)


@pytest.fixture(scope="module")
def reference():
    return load(IG_PDB)


def test_reference_atom_counts(reference):
    assert reference.chains == ["A", "B", "C", "D"]
    for chain, n in ATOMS_PER_CHAIN.items():
        assert len(load(IG_PDB, chain=chain).coord) == n
    assert len(reference.coord) == sum(ATOMS_PER_CHAIN.values())


@pytest.mark.parametrize("fmt", ["pdb.gz", "cif", "CIF", "mmcif", "cif.gz"])
def test_formats_match_pdb(coord_files, reference, fmt):
    """Every accepted extension (any case, optionally gzipped) loads the same atoms."""
    s = load(coord_files[fmt])
    assert s.chains == reference.chains
    assert s.sequence == reference.sequence
    np.testing.assert_allclose(s.coord, reference.coord, atol=1e-3)


def test_unknown_extension_exits(tmp_path):
    bad = tmp_path / "1igt.xyz"
    bad.write_text("")
    with pytest.raises(SystemExit, match="File format not recognized"):
        load(bad)


@pytest.mark.parametrize("chain, expected", [
    ("all", ["A", "B", "C", "D"]),
    (["all"], ["A", "B", "C", "D"]),
    ("A", ["A"]),
    (["A", "C"], ["A", "C"]),
    ("AC", ["A", "C"]),        # older shorthand: one string of single-letter IDs
    (["AC"], ["A", "C"]),      # what the CLI passes for --chain AC
    (["C", "A", "C"], ["C", "A"]),
])
def test_chain_selection(chain, expected):
    s = load(IG_PDB, chain=chain)
    assert s.chains == expected
    assert len(s.coord) == sum(ATOMS_PER_CHAIN[c] for c in expected)


def test_unknown_chain_exits():
    with pytest.raises(SystemExit, match="Chain 'Z' not found"):
        load(IG_PDB, chain=["Z"])


def test_multi_letter_chain_ids(coord_files):
    path = coord_files["multi"]
    assert load(path).chains == ["AA", "BB", "C", "D"]

    # "AA" is matched as one chain, not split into A + A
    s = load(path, chain=["AA"])
    assert s.chains == ["AA"]
    assert len(s.coord) == ATOMS_PER_CHAIN["A"]

    s = load(path, chain=["BB", "C"])
    assert s.chains == ["BB", "C"]
    assert len(s.coord) == ATOMS_PER_CHAIN["B"] + ATOMS_PER_CHAIN["C"]

    # single-letter shorthand still works alongside multi-letter IDs
    assert load(path, chain="CD").chains == ["C", "D"]


def test_multi_letter_chain_ids_match_renamed_pdb_chains(coord_files, reference):
    s = load(coord_files["multi"])
    for old, new in MULTI_LETTER_CHAINS.items():
        assert s.sequence[new] == reference.sequence[old]
