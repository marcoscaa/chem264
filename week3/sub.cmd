#!/bin/bash

#SBATCH -p instruction  # Partition name
#SBATCH -J h2o_md        # Job name
#SBATCH --mail-user=<cruzid>@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH -o job%.j.out    # Name of stdout output file
#SBATCH -n 4        # Total number of mpi tasks requested per node
#SBATCH -t 01:00:00  # Run Time (hh:mm:ss) - 30 min (optional)
#SBATCH --mem=1G # Memory to be allocated PER NODE

export OMPI_MCA_btl=tcp,sm,self
module load quantumespresso/7.2

# Use of -p replaces the need to use "#SBATCH --cpus-per-task"
mpirun -np $SLURM_NTASKS pw.x < h2o_nvt.in > h2o_nvt.out
