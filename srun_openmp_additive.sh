#!/usr/bin/bash
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 12:00:00
#SBATCH -A mor100
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH --export=ALL

# Parallel programs run on a dedicated node using all (or at least half) of the cores

NUM_THREADS=128  # 64 may be a good choice as well

SCRATCH="/scratch/$USER/job_$SLURM_JOB_ID"

cp "$HOME/hpc-final/matrix.c" "$SCRATCH"
cp "$HOME/hpc-final/matrix.h" "$SCRATCH"
cp "$HOME/hpc-final/matrix_io_helpers.h" "$SCRATCH"
cp "$HOME/hpc-final/util.c" "$SCRATCH"
cp "$HOME/hpc-final/util.h" "$SCRATCH"

FILE="2d-fdtd-additive-p.c"
cp "$HOME/hpc-final/$FILE" "$SCRATCH"
SFILE="$SCRATCH/$FILE"

echo "additive serial"

for n in 100 500 1000; do # 5000 10000 500000; do
    echo $n
    gcc -Wall -g -O3 -march=native $SFILE matrix.c util.c -o $SCRATCH/2d-fdtd-additive-p -lm -DMESH_SIZE=$n -DNUM_TIMESTEPS=500 -DSAVE_EVERY_N_STEPS=5
    FNAME="output-fdtd-2d-additive-p-$n.npy"
    "$SCRATCH/2d-fdtd-additive-p" "$FNAME" # & means to background the process
    "$SCRATCH/2d-fdtd-additive-p" "$FNAME"
    "$SCRATCH/2d-fdtd-additive-p" "$FNAME"
done

