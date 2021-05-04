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

for PROGNAME in transverse additive interval; do
    echo 2d-fdtd-$PROGNAME-s
    cp "$HOME/hpc-final/2d-fdtd-$PROGNAME-s.c" "$SCRATCH"
    SFILE="$SCRATCH/2d-fdtd-$PROGNAME-s.c"
    TIMESTEPS=100
    for MESH_SIZE in {100..1000..$TIMESTEPS} # allowed?
    # for MESH_SIZE in {50..100..10}
    do
        echo "SIZE: $MESH_SIZE"
        echo "TIMESTEPS: $TIMESTEPS"
        gcc -Wall -g -O3 -march=native $SFILE matrix.c util.c -o $SCRATCH/2d-fdtd-$PROGNAME-s -lm -DMESH_SIZE=$MESH_SIZE -DNUM_TIMESTEPS=$TIMESTEPS -DSAVE_EVERY_N_STEPS=2
        FNAME="output-fdtd-2d-$PROGNAME-$MESH_SIZE-s.npy"
        "$SCRATCH/2d-fdtd-$PROGNAME-s" "$FNAME" # & means to background the process
        "$SCRATCH/2d-fdtd-$PROGNAME-s" "$FNAME"
        "$SCRATCH/2d-fdtd-$PROGNAME-s" "$FNAME"
        "$SCRATCH/2d-fdtd-$PROGNAME-s" "$FNAME"
        "$SCRATCH/2d-fdtd-$PROGNAME-s" "$FNAME"
    done
    #wait # waits for all background processes
done

