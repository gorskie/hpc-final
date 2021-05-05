/**
 * To compile:  gcc -Wall -g -O3 -fopenmp 2d-fdtd-additive-p.c matrix.c util.c -o 2d-fdtd-additive-p -lm
 * To run: ./2d-fdtd-additive-p
**/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "matrix.h"
#include "util.h"

#ifndef MESH_SIZE
#define MESH_SIZE 2000
#endif
#define MESH_SIZE_SQUARED MESH_SIZE*MESH_SIZE
#ifndef NUM_TIMESTEPS
#define NUM_TIMESTEPS 500
#endif
#ifndef SAVE_EVERY_N_STEPS
#define SAVE_EVERY_N_STEPS 10
#endif
#ifndef T0
#define T0 20 // >= 2
#endif
#ifndef PULSE_SPREAD
#define PULSE_SPREAD 6.f // gaussian
#endif
#ifndef PULSE_WIDTH
#define PULSE_WIDTH 2*MESH_SIZE/3
#endif
#ifndef PULSE_Y
#define PULSE_Y 2 // apply slightly below the boundary
#endif
#ifndef PULSE_X
#define PULSE_X MESH_SIZE/2 // apply in middle of the mesh
#endif
#ifndef CELL_SIZE
#define CELL_SIZE 3e8f
#endif
#ifndef DZ
#define DZ  0.01f
#endif
#ifndef DT
#define DT  DZ/(2*CELL_SIZE)
#endif
#ifndef CC
#define CC CELL_SIZE*DT/DZ // 0.5f unless DT is not the derived value that it is
#endif

#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define MAX(X,Y) ((X)>(Y)?(X):(Y))

int main(int argc, const char *argv[]) {
    // parse args
    const char* output_loc = "./output-fdtd-2d-additive-p.npy";
    size_t num_threads = 4;

    if(argc > 1) {
        num_threads = argv[1];
    }
    if (argc > 2) {
        output_loc = argv[2];
    }


    float* ez = (float*)malloc(MESH_SIZE_SQUARED*sizeof(float));
    float* hx = (float*)malloc(MESH_SIZE_SQUARED*sizeof(float));
    float* hy = (float*)malloc(MESH_SIZE_SQUARED*sizeof(float));
    float *ez_row, *hy_row, *hx_row;
    Matrix *output = matrix_zeros(
        (NUM_TIMESTEPS+SAVE_EVERY_N_STEPS-1)/SAVE_EVERY_N_STEPS,
        MESH_SIZE*MESH_SIZE);

    // start the clock
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // CHECK: this correct?
    #pragma omp parallel num_thread(num_threads)
    for (int t = 0, steps_until_printout = SAVE_EVERY_N_STEPS; t < NUM_TIMESTEPS; t++) {

        /* --- MAIN FDTD LOOP --- */
        // TODO: check if need loop fusion
        ez_row = ez; hy_row = hy; hx_row = hx;
        #pragma omp for
        for (int i = 1; i < MESH_SIZE; i++) {
            ez_row += MESH_SIZE; hy_row += MESH_SIZE; hx_row += MESH_SIZE;
            for (int j = 1; j < MESH_SIZE; j++) {
                ez_row[j] += CC*(hy_row[j] - hy_row[j-MESH_SIZE] - hx_row[j] + hx_row[j-1]);
            }
        }

        // pulse centered at PULSE_X, PULSE_Y with PULSE_WIDTH width and PULSE_SPREAD "height"
        // adds pulse every timestep at a constant intensity
        float pulse = 1.f / PULSE_SPREAD;
        pulse = exp(-0.5f * pulse * pulse) * (2.f / PULSE_WIDTH);
        ez_row = ez+PULSE_Y*MESH_SIZE;
        #pragma omp for
        for (int j = MAX(PULSE_X-PULSE_WIDTH/2, 0); j < MIN(PULSE_X+PULSE_WIDTH/2, MESH_SIZE); j++) {
            ez_row[j] += pulse;
        }
        
        // TODO: loop fission
        // Calculate magnetic field in the X
        // Calculate magnetic field in the Y
        ez_row = ez; hy_row = hy; hx_row = hx;
        #pragma omp for
        for (int i = 1; i < MESH_SIZE-1; i++) {
            for (int j = 1; j < MESH_SIZE-1; j++) {
                hx_row[j] += CC*(ez_row[j] - ez_row[j+1]);
                hy_row[j] += CC*(ez_row[MESH_SIZE+j] - ez_row[j]);
            }
            ez_row += MESH_SIZE; hy_row += MESH_SIZE; hx_row += MESH_SIZE;
        }

        /* --- END OF MAIN FDTD LOOP --- */

        if (--steps_until_printout == 0) {
            memcpy(output->data + t/SAVE_EVERY_N_STEPS*MESH_SIZE_SQUARED, ez, MESH_SIZE_SQUARED*sizeof(float));
            steps_until_printout = SAVE_EVERY_N_STEPS;
        }
    }
    // get the end and computation time
    clock_gettime(CLOCK_MONOTONIC, &end);
    double time = get_time_diff(&start, &end);
    printf("%f\n", time);

    // save results
    matrix_to_npy_path(output_loc, output);
    // cleanup
    free(ez);
    free(hx);
    free(hy);
    matrix_free(output);
    return 0;
}
