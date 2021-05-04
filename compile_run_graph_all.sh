#!/bin/bash
gcc -Wall -g -O3 -march=native 2d-fdtd-transverse-s.c matrix.c util.c -o 2d-fdtd-transverse-s -lm -DMESH_SIZE=100 -DNUM_TIMESTEPS=251 -DSAVE_EVERY_N_STEPS=5
gcc -Wall -g -O3 -march=native 2d-fdtd-additive-s.c matrix.c util.c -o 2d-fdtd-additive-s -lm -DMESH_SIZE=100 -DNUM_TIMESTEPS=251 -DSAVE_EVERY_N_STEPS=5
gcc -Wall -g -O3 -march=native 2d-fdtd-interval-s.c matrix.c util.c -o 2d-fdtd-interval-s -lm -DMESH_SIZE=100 -DNUM_TIMESTEPS=251 -DSAVE_EVERY_N_STEPS=5
./2d-fdtd-transverse-s
./2d-fdtd-additive-s
./2d-fdtd-interval-s
python print-graph.py 100 transverse
python print-graph.py 100 additive
python print-graph.py 100 interval