for file in scf_*.out 
do 
  alat=$(echo "$file" | sed 's/scf_\([0-9.]*\)\.out/\1/')
  a=`grep ! scf_${alat}.out | awk '{print $5}'` 
  echo $alat $a 
done > data

cat << EOF > plot.py
import numpy as np
import matplotlib.pyplot as plt
f=np.loadtxt("data").T
plt.figure(figsize=(3,3))
plt.plot(f[0],f[1])
plt.xlabel("alat (\$\\AA\$)")
plt.ylabel("Energy (Ry)")
plt.tight_layout()
plt.savefig('cell_relax.pdf')
EOF

python plot.py
rm plot.py data
