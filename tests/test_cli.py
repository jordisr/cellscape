"""Smoke tests that drive the ``cellscape`` CLI end to end on the bundled examples."""
import pickle
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

from cellscape import cli
from conftest import IG_PDB, IG_VIEW


def run_cli(monkeypatch, *argv):
    monkeypatch.setattr(sys, "argv", ["cellscape", *map(str, argv)])
    try:
        cli.main()
    finally:
        plt.close("all")


def export_cartoon(monkeypatch, out_dir, name, *extra, coords=IG_PDB):
    """Export a cartoon pickle (alpha carbons only, to keep the tests fast)."""
    stem = out_dir / name
    run_cli(monkeypatch, "cartoon", "--pdb", coords, "--view", IG_VIEW, "--only_ca",
            "--export", "--save", f"{stem}.png", *extra)
    pickle_path = Path(f"{stem}.pickle")
    assert pickle_path.exists()
    return pickle_path


@pytest.fixture(scope="module")
def cartoons(tmp_path_factory):
    """Cartoon pickles with and without depth shading, built once per module."""
    out_dir = tmp_path_factory.mktemp("cartoons")
    with pytest.MonkeyPatch.context() as mp:
        return {
            "flat": export_cartoon(mp, out_dir, "flat", "--outline", "chain", "--depth", "flat"),
            "shaded": export_cartoon(mp, out_dir, "shaded", "--outline", "residue",
                                     "--color_by", "chain", "--depth_shading"),
        }


@pytest.mark.parametrize("extra", [
    ["--outline", "all", "--back_outline"],
    ["--outline", "chain", "--depth", "flat", "--back_outline"],
    ["--outline", "chain", "--depth", "contours", "--back_outline"],
    ["--outline", "residue", "--color_by", "chain", "--depth_shading", "--depth_lines"],
])
def test_cartoon_renders(monkeypatch, tmp_path, extra):
    export_cartoon(monkeypatch, tmp_path, "cartoon", *extra)
    assert (tmp_path / "cartoon.png").stat().st_size > 0


def test_scene_renders(monkeypatch, tmp_path, cartoons):
    out = tmp_path / "scene.png"
    run_cli(monkeypatch, "scene", "--files", cartoons["flat"], cartoons["shaded"], "--save", out)
    assert out.stat().st_size > 0


@pytest.mark.parametrize("which", ["flat", "shaded"])
def test_scene_recolor(monkeypatch, tmp_path, cartoons, which):
    # Cartoons exported without --depth_shading store "shade": None for each
    # polygon; recoloring used to crash on them with a TypeError.
    out = tmp_path / "scene.png"
    run_cli(monkeypatch, "scene", "--files", cartoons[which], "--recolor",
            "--recolor_cmap", "viridis", "--save", out)
    assert out.stat().st_size > 0


def chain_areas(pickle_path):
    with open(pickle_path, "rb") as f:
        return [round(p["polygon"].area, 1) for p in pickle.load(f)["polygons"]]


@pytest.mark.parametrize("fmt", ["cif", "cif.gz"])
def test_cartoon_from_cif_matches_pdb(monkeypatch, tmp_path, coord_files, fmt):
    args = ["--outline", "chain", "--depth", "flat"]
    from_pdb = export_cartoon(monkeypatch, tmp_path, "from_pdb", *args)
    from_cif = export_cartoon(monkeypatch, tmp_path, "from_cif", *args, coords=coord_files[fmt])
    assert chain_areas(from_cif) == chain_areas(from_pdb)


def test_cartoon_multi_letter_chain(monkeypatch, tmp_path, coord_files):
    # --chain AA must select the one chain named "AA", matching chain A of the PDB
    from_pdb = export_cartoon(monkeypatch, tmp_path, "chain_a", "--chain", "A")
    from_cif = export_cartoon(monkeypatch, tmp_path, "chain_aa", "--chain", "AA",
                              coords=coord_files["multi"])
    assert chain_areas(from_cif) == chain_areas(from_pdb)
