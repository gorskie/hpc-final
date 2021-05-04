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

for PROGNAME in transverse additive interval; do
    echo $PROGNAME
    cp "$HOME/2d-fdtd-$PROGNAME-p.c" "$SCRATCH"
    SFILE="$SCRATCH/2d-fdtd-$PROGNAME-p.c"
    for MESH_SIZE in {50..100..10}
    do  
        echo $MESH_SIZE
        gcc -Wall -g -O3 -march=native 2d-fdtd-$PROGNAME-p.c matrix.c util.c -o 2d-fdtd-$PROGNAME-p -lm -DMESH_SIZE=$MESH_SIZE
        FNAME="output-fdtd-2d-$FILE-$MESH_SIZE-p.npy"
        ./2d-fdtd-$PROGNAME-p "$FNAME" &  # & means to background the process
    done
    wait # waits for all background processes
    cp "$SCRATCH/*.npy" "$HOME/"
done