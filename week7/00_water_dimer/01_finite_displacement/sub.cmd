#!/bin/bash
#SBATCH -p instruction  # Partition name
#SBATCH -J test        # Job name
#SBATCH --mail-type=ALL
#SBATCH -o job%.j.out    # Name of stdout output file
#SBATCH -N 1        # Total number of nodes requested (128x24/Instructional only)
#SBATCH -n 16        # Total number of mpi tasks requested per node
#SBATCH -t 00:30:00  # Run Time (hh:mm:ss) - 30 min (optional)
#SBATCH --mem=2G # Memory to be allocated PER NODE

export OMPI_MCA_btl=tcp,sm,self
module load quantumespresso/7.2

# I am assuming that you generated the input files using the instructions in README.md. If you have not done so, please do that first. 

for i in *.in
do
    echo "Running $i"
    mpirun -np $SLURM_NTASKS pw.x -nk 4 < $i > ${i%.in}.out
done