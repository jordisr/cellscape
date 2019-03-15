# Tools for proteome visualization
## pdb2svg: Vector outlines of macromolecular structure
<img src="ig_example.png" alt="logo" width=700/>

### Requirements
To run `pdb2svg.py` you will need...
* Python 3 (also biopython and shapely libraries)
* [PyMOL](https://pymol.org/2/) (optional)

If you have Python 3 installed you should be able to get the dependencies with
```
pip install biopython
pip install shapely
```

### Example: 1IGT
We can download an immunoglobulin structure from the PDB to test on (used to generate the above images):
```
curl -O https://files.rcsb.org/view/1IGT.pdb
```
#### Selecting the camera view in PyMOL
First open the protein structure in PyMOL, choose the desired rotation (zoom is irrelevant), and enter `get_view` in the PyMOL console. The output should look something like this:
```
### cut below here and paste into script ###
set_view (\
     0.761977673,    0.134355009,    0.633512199,\
    -0.555662155,    0.638066173,    0.533020258,\
    -0.332609177,   -0.758168101,    0.560847342,\
    -0.000003876,    0.000002533, -7417.434570312,\
     0.751698136,  -16.697633743,    9.022263527,\
  7324.144042969, 7506.925292969,  -20.000000000 )
### cut above here and paste into script ###
```
Copy and paste the indicated region into a new text file, e.g. `view.txt`.

#### Generating graphics
...
