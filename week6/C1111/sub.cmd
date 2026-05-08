#!/bin/bash
#SBATCH -p instruction  # Partition name
#SBATCH -J test        # Job name
#SBATCH --mail-type=ALL
#SBATCH -o job%.j.out    # Name of stdout output file
#SBATCH -N 1        # Total number of nodes requested (128x24/Instructional only)
#SBATCH -n 16        # Total number of mpi tasks requested per node
#SBATCH -t 00:30:00  # Run Time (hh:mm:ss) - 30 min (optional)
#SBATCH --mem=4G # Memory to be allocated PER NODE

export OMPI_MCA_btl=tcp,sm,self
module load quantumespresso/7.2

cp ../../Ge.upf .

mpirun -np $SLURM_NTASKS pw.x -nk 4 < plus_e.in  > plus_e.out
mpirun -np $SLURM_NTASKS pw.x -nk 4 < minus_e.in > minus_e.out
