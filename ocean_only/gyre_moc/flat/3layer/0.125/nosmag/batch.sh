#!/bin/bash
#
#SBATCH --job-name=flat_0.125
#SBATCH --output=output.txt
#SBATCH --ntasks-per-node=28
#SBATCH --nodes=4
#SBATCH --time=48:00:00
#SBATCH -p long-28core
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=christopher.wolfe@stonybrook.edu

module load intel/compiler/64/2020/20.0.2
module load intel/mkl/64/2020/20.0.2
module load intel/mpi/64/2020/20.0.2
module load netcdf-fortran

cd $SLURM_SUBMIT_DIR
mkdir -p OUTPUT
mkdir -p RESTART

date
mpirun ./MOM6
date

