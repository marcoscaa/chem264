for sys in 01_change_alat 02_change_alat_beta_tin 
do
  cd $sys
  for file in scf_*.out 
  do 
    alat=$(echo "$file" | sed 's/scf_\([0-9.]*\)\.out/\1/')
    a=`grep ! scf_${alat}.out | awk '{print $5}'` 
    v=`grep volume scf_${alat}.out | awk '{print $4}'`
    echo $alat $a $v 
  done > ../data_$sys
  cd ..
done

cat << EOF > plot.py
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

p=float(argv[1])
e_iso=-178.37283191
ry2ev=13.6057
au2GPa=ry2ev/0.529177**3*160.2
f=np.loadtxt("data_01_change_alat").T
g=np.loadtxt("data_02_change_alat_beta_tin").T
f[1]=(f[1]-2*e_iso)/2.
g[1]=(g[1]-2*e_iso)/2.

plt.figure(figsize=(3,3))
plt.plot(f[0],f[1]*ry2ev,label="Diamond")
plt.plot(g[0],g[1]*ry2ev,label="beta-Tin")
plt.xlabel("alat (\$\\AA\$)")
plt.ylabel("Cohesive energy (eV)")
plt.legend()
plt.tight_layout()
plt.savefig('cell_relax.pdf')

h1=(f[1]+p*f[2]/au2GPa/2.)
h2=(g[1]+p*g[2]/au2GPa/2.)
plt.figure(figsize=(3,3))
plt.plot(f[0],h1*ry2ev,label="Diamond")
plt.plot(g[0],h2*ry2ev,label="beta-Tin")
plt.xlabel("alat (\$\\AA\$)")
plt.ylabel("Cohesive enthalpy (eV)")
plt.legend()
plt.tight_layout()
plt.savefig('cell_relax_enthalpy.pdf')
EOF

python plot.py 10
rm plot.py data_01_change_alat data_02_change_alat_beta_tin
