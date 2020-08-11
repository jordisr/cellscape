# CellScape: Protein structure visualization with vector graphics cartoons
<img src="ig_example.png" alt="logo" width=700/>

## Installation
To run CellScape you will need:
* Python 3
* [PyMOL](https://pymol.org/2/) or [Chimera](https://www.cgl.ucsf.edu/chimera/) (optional, needed to orient the protein if not using the Jupyter notebook interface)

CellScape and its dependencies can be installed with:

```
git clone https://github.com/jordisr/cellscape
cd cellscape
pip install -r requirements.txt
pip install -e .
```

## Making a cartoon from a PDB structure (`cellscape cartoon`)

### Jupyter notebook interface
The most interactive way of building cartoons is through the Python package interface. An example notebook is provided [here](examples/cartoon.ipynb).

### Command-line interface

Cartoons can also be built in one-go from the command-line, as illustrated below.

#### Generating molecular outlines
The following examples should yield the three images used in the top figure (from left to right):
```
cellscape cartoon --pdb examples/ig/1igt.pdb --view examples/ig/view --outline residue --color_by chain
```
The most realistic visualization projects the 3D coordinates down to two dimensions, and outlines each residue separately. Shading is used to simulate depth in a style inspired by [David Goodsell](https://pdb101.rcsb.org/motm/21).

```
cellscape cartoon --pdb examples/ig/1igt.pdb --view examples/ig/view --outline chain --depth flat
```
Each chain is outlined separately. The `--depth flat` option ensures that if the chains overlap, only the portion that is visible (i.e. closer to the camera) is incorporated into the outline.

```
cellscape cartoon --pdb examples/ig/1igt.pdb --view examples/ig/view --outline all
```
A simple space-filling outline of the entire structure.

Full description of all options is available by running `cellscape cartoon -h`.

### Exporting the camera view
The camera orientation can be set interactively through the Jupyter notebook interface, however to use the command-line interface you will need a separate file with the rotation matrix.
One option is to export it from another molecular visualization tool (currently PyMOL and Chimera formats are supported).

#### PyMOL
Open the protein structure in PyMOL, and choose the desired rotation (zoom is irrelevant). Next, enter `get_view` in the PyMOL console. The output should look something like this:
```
### cut below here and paste into script ###
set_view (\
    -0.273240060,   -0.516133010,    0.811750829,\
     0.870557129,    0.226309016,    0.436930388,\
    -0.409222305,    0.826064587,    0.387488008,\
     0.000000000,    0.000000000, -544.673034668,\
    -0.071666718,  -17.390396118,    8.293336868,\
   455.182373047,  634.163574219,  -20.000000000 )
### cut above here and paste into script ###
```
Copy and paste the indicated region (between the ### lines) into a new text file, which can be passed to CellScape.

#### Chimera
Open the protein structure in Chimera, and choose the desired rotation (zoom is irrelevant).
Enter the command `matrixget` (if no output filename is given it will prompt you for one).
This will write the rotation matrix to a file that can be understood by CellScape.
It should look something like this:
```
Model 0.0
        -0.607365 0.792409 0.0565265 9.04218
        -0.309318 -0.301425 0.901923 -30.7393
        0.731731 0.530312 0.428181 15.789
```

## Composing cartoons into a cellular scene (`cellscape scene`)

Full description of all options is available by running `cellscape scene -h`.
