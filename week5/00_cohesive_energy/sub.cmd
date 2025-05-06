#!/bin/bash

#SBATCH -p instruction  # Partition name
#SBATCH -J test        # Job name
#SBATCH --mail-user=<cruzid>@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH -o job%.j.out    # Name of stdout output file
#SBATCH -N 1        # Total number of nodes requested (128x24/Instructional only)
#SBATCH -n 16        # Total number of mpi tasks requested per node
#SBATCH -t 00:30:00  # Run Time (hh:mm:ss) - 30 min (optional)
#SBATCH --mem=2G # Memory to be allocated PER NODE

export OMPI_MCA_btl=tcp,sm,self
module load quantumespresso/7.2

# Use of -p replaces the need to use "#SBATCH --cpus-per-task"
mpirun -np $SLURM_NTASKS pw.x -nk 4 < vc-relax_bulk_diamond.in > vc-relax_bulk_diamond.out 
mpirun -np $SLURM_NTASKS pw.x -nk 4 < vc-relax_bulk_b-tin.in > vc-relax_bulk_b-tin.out
mpirun -np 4 pw.x < ge_isolated.in > ge_isolated.out
