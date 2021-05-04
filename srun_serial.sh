#!/usr/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=10G
#SBATCH -t 12:00:00
#SBATCH -A mor100
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH --export=ALL

# Serial programs run on shared nodes using a few cores to run individual processes in parallel

SCRATCH="/scratch/$USER/job_$SLURM_JOB_ID"

for PROGNAME in transverse additive interval; do
    echo $PROGNAME
    cp "$HOME/2d-fdtd-$PROGNAME-s.c" "$SCRATCH"
    SFILE="$SCRATCH/2d-fdtd-$PROGNAME-s.c"
    for MESH_SIZE in {50..100..10}
    do  
        echo $MESH_SIZE
        gcc -Wall -g -O3 -march=native 2d-fdtd-$PROGNAME-s.c matrix.c util.c -o 2d-fdtd-$PROGNAME-s -lm -DMESH_SIZE=$MESH_SIZE
        FNAME="output-fdtd-2d-$FILE-$MESH_SIZE-s.npy"
        ./2d-fdtd-$PROGNAME-s "$FNAME" &  # & means to background the process
    done
    wait # waits for all background processes
    cp "$SCRATCH/*.npy" "$HOME/"
done
