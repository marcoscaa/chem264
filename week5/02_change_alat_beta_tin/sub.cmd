#!/bin/bash

#SBATCH -p instruction  # Partition name
#SBATCH -J test        # Job name
#SBATCH --mail-user=<cruzid>@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH -o job%.j.out    # Name of stdout output file
#SBATCH -N 1        # Total number of nodes requested (128x24/Instructional only)
#SBATCH -n 16        # Total number of mpi tasks requested per node
#SBATCH -t 02:00:00  # Run Time (hh:mm:ss) - 30 min (optional)
#SBATCH --mem=2G # Memory to be allocated PER NODE

export OMPI_MCA_btl=tcp,sm,self
module load quantumespresso/7.2

# Use of -p replaces the need to use "#SBATCH --cpus-per-task"

for alat in "4.0" "4.2" "4.4" "4.6" "4.8" "5.0" "5.1" "5.2" "5.3" "5.4" "5.6" "5.8" 
do
  #celldm=`echo "${alat}*1.88973" | bc -l`
  celldm=`python -c "print(${alat}*1.88973)"`
  echo $celldm
  sed "s/CHANGEME/${celldm}/g" scf_t.in > scf.in 
  mpirun -np $SLURM_NTASKS pw.x -nk 4 < scf.in > scf_${alat}.out
done
