import argparse, sys, os
from .cartoon import make_cartoon
from .scene import make_scene

def main():
    # set up argument parser
    parser = argparse.ArgumentParser(description='CellScape: Vector graphics for molecular visualization')
    subparsers = parser.add_subparsers(dest="command")
    subparsers.required=True

    # cartoon
    parser_cartoon = subparsers.add_parser('cartoon', help='', formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="")
    parser_cartoon.set_defaults(func=make_cartoon)
    # input/output options
    parser_cartoon_io = parser_cartoon.add_argument_group('input/output options')
    parser_cartoon_io.add_argument('--pdb', help='Protein coordinates file (must be .pdb/.ent/.cif/.mcif)', required=True)
    parser_cartoon_io.add_argument('--model', type=int, default=0, help='Model number in PDB to load')
    parser_cartoon_io.add_argument('--chain', default=['all'], help='Chain(s) in structure to outline', nargs='+')
    parser_cartoon_io.add_argument('--view', help='Camera rotation matrix (saved from cellscape, PyMOL get_view, or Chimera matrixget)')
    parser_cartoon_io.add_argument('--uniprot', help='UniProt XML file to parse for sequence/domain/topology information')
    parser_cartoon_io.add_argument('--save', default='out.svg', help='Image output file (valid formats are png/pdf/svg/ps)')
    parser_cartoon_io.add_argument('--export', default=False, action="store_true", help='Export Python object with structural information')
    # outline building options
    parser_cartoon_outline = parser_cartoon.add_argument_group('outline-building options')
    parser_cartoon_outline.add_argument('--only_annotated', action='store_true', default=False, help='Ignore regions without UniProt annotations')
    parser_cartoon_outline.add_argument('--only_ca', action='store_true', default=False, help='Only use alpha carbons for outline')
    parser_cartoon_outline.add_argument('--outline_by', '--outline',  default='all',  choices=['all', 'chain', 'domain', 'topology', 'residue'], help='Outline protein regions')
    parser_cartoon_outline.add_argument('--depth',  default=None,  choices=['flat', 'contours', None], help='Represent depth with flat occluded outlines or contour slices')
    parser_cartoon_outline.add_argument('--depth_contour_interval',  default=3, help='Width of depth contour bins in angstroms (if --depth contours)')
    parser_cartoon_outline.add_argument('--radius', default=1.5, help='Atomic radius, in angstroms', type=float)
    # visual style options
    parser_cartoon_style = parser_cartoon.add_argument_group('styling options')
    parser_cartoon_style.add_argument('--axes', action='store_true', default=False, help='Draw x and y axes around molecule')
    parser_cartoon_style.add_argument('--colors', default=[], nargs='+', help='Specify color scheme for protein (list of colors or matplotlib named color map)')
    parser_cartoon_style.add_argument('--edge_color', default='black', help='Edge color')
    parser_cartoon_style.add_argument('--line_width', default=0.7, type=float, help='Line width')
    parser_cartoon_style.add_argument('--color_by', default='same',  choices=['same', 'chain', 'domain', 'topology'], help='Color residues by attribute (if --outline_by residues is selected)')
    parser_cartoon_style.add_argument('--depth_shading', action='store_true', default=False, help='Shade regions darker in the back to simulate depth')
    parser_cartoon_style.add_argument('--depth_lines', action='store_true', default=False, help='Use thicker lines the back to simulate depth')
    parser_cartoon_style.add_argument('--dpi', type=int, default=300, help='DPI to use if exporting to a raster format like PNG')

    # scene
    parser_scene = subparsers.add_parser('scene', help='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_scene.set_defaults(func=make_scene)
    # input/output options
    parser_scene_io = parser_scene.add_argument_group('input/output options')
    parser_scene_io.add_argument('--files', nargs='+', help='Pickled objects to load')
    parser_scene_io.add_argument('--save', default='out.svg', help='Image output path (valid formats are png/pdf/svg/ps)')
    # visual style options
    parser_scene_style = parser_scene.add_argument_group('styling options')
    parser_scene_style.add_argument('--offsets', nargs='+', default=[], help='Vertical offsets for each molecule specified manually')
    parser_scene_style.add_argument('--padding', type=int, default=0, help='Horizontal padding to add between each molecule (in angstroms)')
    parser_scene_style.add_argument('--axes', action='store_true', default=False, help='Draw x and y axes')
    parser_scene_style.add_argument('--membrane', default=None, choices=[None, 'arc', 'flat', 'wave'], help='Draw membrane on X axis')
    parser_scene_style.add_argument('--membrane_thickness', default=40, type=float, help='Thickness of the membrane (in angstroms)')
    parser_scene_style.add_argument('--membrane_lipids', action='store_true', help='Draw lipid head groups')
    parser_scene_style.add_argument('--no_membrane_offset', action='store_true', help=argparse.SUPPRESS) # don't adjust y-axis to position bottom of structure in membrane
    parser_scene_style.add_argument('--order_by', default='input', choices=['input', 'random', 'height','top'], help='How to order proteins in scene')
    parser_scene_style.add_argument('--recolor', action='store_true', default=False, help='Recolor proteins in scene')
    parser_scene_style.add_argument('--recolor_cmap', default=['hsv'], nargs='+', help='Named cmap or color scheme for re-coloring')
    parser_scene_style.add_argument('--dpi', type=int, default=300, help='DPI to use if exporting to a raster format like PNG')
    parser_scene_style.add_argument('--use_placeholders', action='store_true', help=argparse.SUPPRESS)
    parser_scene_style.add_argument('--labels', action='store_true', default=False, help=argparse.SUPPRESS) # still testing
    parser_scene_style.add_argument('--label_size', type=float, default=0.5, help=argparse.SUPPRESS) # fraction of the screen to use for labels
    parser_scene_style.add_argument('--label_orientation', choices=["vertical", "horizontal", "diagonal"], default="vertical", help=argparse.SUPPRESS)
    parser_scene_style.add_argument('--label_position', choices=["above", "below"], default="below", help=argparse.SUPPRESS)
    parser_scene_style.add_argument('--fig_height', type=float, default=11, help=argparse.SUPPRESS) # passed to figsize
    parser_scene_style.add_argument('--fig_width', type=float, default=8.5, help=argparse.SUPPRESS) # passed to figsize
    # for simulating according to stoichiometry
    parser_scene_sim = parser_scene.add_argument_group('random scene options')
    parser_scene_sim.add_argument('--csv', help='Table of protein information')
    parser_scene_sim.add_argument('--seed', type=int, help='Random seed for scene generation')
    parser_scene_sim.add_argument('--sample_from', help='Column to use for sampling (with --csv)', default='stoichiometry')
    parser_scene_sim.add_argument('--num_mol', type=int, help='Number of molecules to sample for scene', default=0)
    parser_scene_sim.add_argument('--background', action='store_true', default=False, help='Add background plane using same frequencies')

    # parse arguments and call corresponding command
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
