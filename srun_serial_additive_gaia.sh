#!/usr/bin/bash
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH --export=ALL

# Serial programs run on shared nodes using a few cores to run individual processes in parallel

SCRATCH="/scratch/$USER/job_$SLURM_JOB_ID"

cp "$HOME/hpc-final/matrix.c" "$SCRATCH"
cp "$HOME/hpc-final/matrix.h" "$SCRATCH"
cp "$HOME/hpc-final/matrix_io_helpers.h" "$SCRATCH"
cp "$HOME/hpc-final/util.c" "$SCRATCH"
cp "$HOME/hpc-final/util.h" "$SCRATCH"

FILE="2d-fdtd-additive-s.c"
cp "$HOME/hpc-final/$FILE" "$SCRATCH"
SFILE="$SCRATCH/$FILE"

echo "additive serial"

for n in 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000; do
    echo $n
    gcc -Wall -g -O3 -march=native $SFILE matrix.c util.c -o $SCRATCH/2d-fdtd-additive-s -lm -DMESH_SIZE=$n -DNUM_TIMESTEPS=100 -DSAVE_EVERY_N_STEPS=5
    FNAME="output-fdtd-2d-additive-s-$n.npy"
    "$SCRATCH/2d-fdtd-additive-s" "$FNAME" # & means to background the process
    "$SCRATCH/2d-fdtd-additive-s" "$FNAME"
    "$SCRATCH/2d-fdtd-additive-s" "$FNAME"
done

