for alat in "5.2" "5.3" "5.4" "5.5" "5.6" "5.7" "5.72" "5.74" "5.76" "5.78" "5.8" "5.82" "5.84" "5.9" "6.0" "6.1" "6.2" "6.3" "6.4" "6.5" "6.6" 
do 
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
