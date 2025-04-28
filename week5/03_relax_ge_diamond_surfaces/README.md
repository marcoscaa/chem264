# Make surface supercells

You can generate a surface supercell using the code ```create_pw_input_surface_from_materials_project.py``` in this folder. 

## Instructions

```
python create_pw_input_surface_from_materials_project.py Material_ID miller_index nA nB
```

Material_ID: material ID from Materials Project
miller_index: sequence of 3 integers (no spaces). Ex: type 100 for the (100) surface.
nA: replicate the cell nA times along the x direction
nB: replicate the cell nB times along the y direction

This script will generate a PW input file called ```relax_surface.in```. You need to change the number of K-points in this script.
