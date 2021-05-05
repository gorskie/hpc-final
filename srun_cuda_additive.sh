#!/usr/bin/bash
#SBATCH -p gpu-shared
#SBATCH --gpus=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=95G
#SBATCH -t 01:00:00
#SBATCH -A mor100
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH --export=ALL

module purge
module load gpu
module load slurm
module load gcc/9.2.0
module load cuda/11.0.2

SCRATCH="/scratch/$USER/job_$SLURM_JOB_ID"

cp "$HOME/hpc-final/matrix.c" "$SCRATCH"
cp "$HOME/hpc-final/matrix.h" "$SCRATCH"
cp "$HOME/hpc-final/matrix_io_helpers.h" "$SCRATCH"
cp "$HOME/hpc-final/util.c" "$SCRATCH"
cp "$HOME/hpc-final/util.h" "$SCRATCH"

FILE="2d-fdtd-additive.cu"
cp "$HOME/hpc-final/$FILE" "$SCRATCH"
SFILE="$SCRATCH/$FILE"

echo "additive serial"

for n in 100 500 1000 5000 10000 500000; do
    echo $n
    gcc -Wall -O3 -march=native -c matrix.c util.c
    nvcc -arch=sm_20 -O3 $SFILE *.o -o $SCRATCH/2d-fdtd-additive-cuda -lm -DMESH_SIZE=$n -DNUM_TIMESTEPS=500 -DSAVE_EVERY_N_STEPS=5
    FNAME="output-fdtd-2d-additive-cuda-$n.npy"
    "$SCRATCH/2d-fdtd-additive-cuda" "$FNAME"
    "$SCRATCH/2d-fdtd-additive-cuda" "$FNAME"
    "$SCRATCH/2d-fdtd-additive-cuda" "$FNAME"
done

