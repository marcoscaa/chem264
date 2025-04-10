#!/bin/bash

#SBATCH -p instruction  # Partition name
#SBATCH -J test        # Job name
#SBATCH --mail-user=<cruzid>@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH -o job%.j.out    # Name of stdout output file
#SBATCH -N 1        # Total number of nodes requested (128x24/Instructional only)
#SBATCH -n 4        # Total number of mpi tasks requested per node
#SBATCH -t 00:30:00  # Run Time (hh:mm:ss) - 30 min (optional)
#SBATCH --mem=1G # Memory to be allocated PER NODE

export OMPI_MCA_btl=tcp,sm,self
module load quantumespresso/7.2

# Use of -p replaces the need to use "#SBATCH --cpus-per-task"

for alat in "5.2" "5.3" "5.4" "5.5" "5.6" "5.7" "5.72" "5.74" "5.76" "5.78" "5.8" "5.82" "5.84" "5.9" "6.0" "6.1" "6.2" "6.3" "6.4" "6.5" "6.6" 
do
  sed "s/CHAMGEME/${alat}/g" scf_t.in > scf.in 
  mpirun -np $SLURM_NTASKS pw.x < scf.in > scf_${alat}.out
done
