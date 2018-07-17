# CellScape: Tools for proteome visualization
![Logo](cecam5.png)
## pdb2svg
Generate a 2d space-filling outline from a PDB structure.

### To run pdb2svg you will need...
* A PDB structure of your protein of interest
* Python 3 (also biopython and shapely libraries)
* PyMOL

### Overview
1. First open the protein structure in PyMOL, choose the desired rotation (zoom is irrelevant), and enter `get_view` in the PyMOL console. The output should look something like this:
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
Select the indicated region region and paste into a new text file (e.g. `view.txt`).

2. (Optional) If desired, the domain architecture can be specified in a comma-separated file like this:
```
res_start,res_end,description
35,142,Ig-like V-type
145,232,Ig-like C2-type 1
237,319,Ig-like C2-type 2
```

3. Finally run `pdb2svg.py --view view.txt --domains domains.csv` to generate an outline named `out.svg`. Full command-line options are available with `pdb2svg.py --help`.
