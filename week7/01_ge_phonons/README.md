# Computing phonon dispersion using phonopy

First, install phonopy code
```
pip install phonopy
```

Now, create a supercell and apply finite displacements to all atoms in the supercell.
You will need to provide a QE input file and give it to the command below. This 
input file contains the structure of Ge diamond in the conventional unit cell.
```
phonopy --qe -d --dim="2 2 2" --pa=AUTO -c ge_diamond.in
```

This command will generate only 1 PW input file (only the cell and positions). We now have 
to add a header to the file.
```
for file in supercell*in
do
  cat header.in $file > tmp_
  mv tmp_ $file
done
```

Now submit the job
```
sbatch sub.cmd
```

Create the force constants matrix
```
phonopy -f supercell-001.out
```

Post-process the file - plot the phonon dispersion
```
phonopy --qe -c ge_diamond.in -p band.conf
```


