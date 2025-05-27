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
mpirun -np $SLURM_NTASKS pw.x -nk 4 < 00_scf.in > 00_scf.out
mpirun -np $SLURM_NTASKS pw.x -nk 4 < 01_bands.in > 01_bands.out
bands.x < 02_pp_bands.in > 02_pp_bands.out
mpirun -np $SLURM_NTASKS $PW -nk 4 -in 03_dos.in > 03_dos.out
dos.x < 04_pp_dos.in > 04_pp_dos.out

echo ""
echo "Job ended on $(date)"
echo ""

