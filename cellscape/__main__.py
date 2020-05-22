import argparse, sys, os
from .cartoon import make_cartoon
from .scene import make_scene

def main():
    # set up argument parser
    parser = argparse.ArgumentParser(description='CellScape: Vector graphics for molecular visualization')
    subparsers = parser.add_subparsers(dest="command")
    subparsers.required=True

    # cartoon (formerly pdb2svg.py)
    parser_cartoon = subparsers.add_parser('cartoon', help='', formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="WARNING: this software is currently being developed and may not have all the functionality of pdb2svg.py")
    parser_cartoon.set_defaults(func=make_cartoon)
    # input/output options
    parser_cartoon.add_argument('--pdb', help='Input PDB file', required=True)
    parser_cartoon.add_argument('--model', type=int, default=0, help='Model number in PDB to load')
    parser_cartoon.add_argument('--chain', default=['all'], help='Chain(s) in structure to outline', nargs='+')
    parser_cartoon.add_argument('--view', help='File with output from PyMol get_view')
    parser_cartoon.add_argument('--uniprot', nargs='+', help='UniProt XML file to parse for sequence/domain/topology information')
    parser_cartoon.add_argument('--save', default='out', help='Prefix to save graphics')
    parser_cartoon.add_argument('--format', default='svg', help='Format to save graphics', choices=['svg','pdf','png'])
    parser_cartoon.add_argument('--export', action='store_true', help='Export Python object with structural information')
    #parser_cartoon.add_argument('--look', help='Look in directory for structure .pdb, view matrix in .txt and UniProt .xml')
    #parser_cartoon.add_argument('--align', action='store_true', default=False, help='Ignore PDB residue numbering and align to UniProt sequence to find offset')
    parser_cartoon.add_argument('--dpi', type=int, default=300, help='DPI to use if exporting to raster formats (i.e. PNG)')
    parser_cartoon.add_argument('--only_annotated', action='store_true', default=False, help='Ignore regions without UniProt annotations')
    # visual style options
    parser_cartoon.add_argument('--outline_by',  default='all',  choices=['all', 'chain', 'domain', 'topology', 'residue'], help='*')
    parser_cartoon.add_argument('--color_by', default='same',  choices=['same', 'chain', 'domain', 'topology'], help='Color residues by attribute (if --outline_by residues is selected)')
    parser_cartoon.add_argument('--occlude', action='store_true', default=False, help='Occlude residues that are not visible and draw outlines using visible residues only')
    #parser_cartoon.add_argument('--radius', default=1.5, help='Space-filling radius, in angstroms', type=float)
    parser_cartoon.add_argument('--axes', action='store_true', default=False, help='Draw x and y axes around molecule')
    #parser_cartoon.add_argument('--c', default=['#D3D3D3'], nargs='+', help='Set default color(s) in hex RGB')
    #parser_cartoon.add_argument('--cmap', default='Set1', help='Set default color map')
    #parser_cartoon.add_argument('--ec', default='k', help='Set default edge color')
    #parser_cartoon.add_argument('--linewidth', default=0.7, type=float, help='Set default line width')

    # scene (formerly compose_scene.py)
    parser_scene = subparsers.add_parser('scene', help='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_scene.set_defaults(func=make_scene)
    parser_scene.add_argument('--files', nargs='+', help='Pickled objects to load')
    parser_scene.add_argument('--save', default='out', help='Prefix to save graphics')
    parser_scene.add_argument('--offsets', nargs='+', default=[], help='Vertical offsets for each molecule specified manually')
    parser_scene.add_argument('--padding', type=int, default=0, help='Horizontal padding to add between each molecule (in angstroms)')
    parser_scene.add_argument('--axes', action='store_true', default=False, help='Draw x and y axes')
    parser_scene.add_argument('--format', default='png', help='Format to save graphics', choices=['svg','pdf','png'])
    parser_scene.add_argument('--membrane', default=None, choices=[None, 'arc', 'flat', 'wave'], help='Draw membrane on X axis')
    parser_scene.add_argument('--membrane_lipids', action='store_true', help='Draw lipid head groups')
    parser_scene.add_argument('--dpi', type=int, default=300, help='DPI to use if exporting to raster formats (i.e. PNG)')
    parser_scene.add_argument('--order_by', default='input', choices=['input', 'random', 'height'], help='How to order proteins in scene')
    parser_scene.add_argument('--recolor', action='store_true', default=False, help='Recolor proteins in scene')
    parser_scene.add_argument('--recolor_cmap', default=['hsv'], nargs='+', help='Named cmap or color scheme for re-coloring')
    # for simulating according to stoichiometry
    parser_scene.add_argument('--csv', help='Table of protein information')
    parser_scene.add_argument('--sample_from', help='Column to use for sampling', default='stoichiometry')
    parser_scene.add_argument('--num_mol', type=int, help='Total number of molecules in the scene')
    parser_scene.add_argument('--background', action='store_true', default=False, help='Add background plane using same frequencies')

    # parse arguments and call corresponding command
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
