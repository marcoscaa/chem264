# Computing phonon dispersion using phonopy

First, install phono3py code
```
pip install phono3py
```

Now, create a supercell and apply finite displacements to all atoms in the supercell.
You will need to provide a QE input file and give it to the command below. This 
input file contains the structure of Ge diamond in the conventional unit cell.
```
phono3py --qe -d --dim="2 2 2" --pa=AUTO -c ge_diamond.in
```

This command will generate 111 PW input files (only the cell and positions). We now have 
to add a header to each of these files.
```
for file in supercell*.in
do
  cat header.in $file > tmp_$file
  mv tmp_$file $file
done
```

