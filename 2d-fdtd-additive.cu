/**
 * To compile:  gcc -Wall -O3 -march=native -c matrix.c util.c # Run once
 * nvcc -arch=sm_20 -O3 *.cu *.o -o 2d-fdtd-additive-cuda -lm
 * To run: ./2d-fdtd-additive-cuda
**/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

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
#ifndef BLOCK_SIZE
#define BLOCK_SIZE 1024
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

#define CHECK(call)                                                       \
{                                                                         \
   const cudaError_t error = call;                                        \
   if (error != cudaSuccess)                                              \
   {                                                                      \
      printf("Error: %s:%d, ", __FILE__, __LINE__);                       \
      printf("code:%d, reason: %s\n", error, cudaGetErrorString(error));  \
      exit(1);                                                            \
   }                                                                      \
}


__global__
void main_device(float* ez, float* hy, float* hx, float* output) {
    unsigned int i0 = threadIdx.x;
    int steps_until_printout = SAVE_EVERY_N_STEPS;

    for (int t = 0; t < NUM_TIMESTEPS; t++) {
        for (unsigned int i = i0; i < MESH_SIZE*MESH_SIZE; i += blockDim.x) {
            unsigned int x = i % MESH_SIZE, y = i / MESH_SIZE;
            if (x == 0 || y == 0 || x == MESH_SIZE - 1 || y == MESH_SIZE - 1) { continue; }
            float* my_ez = &ez[y*MESH_SIZE + x];
            float* my_hy = &hy[y*MESH_SIZE + x];
            float* my_hx = &hx[y*MESH_SIZE + x];
        
            /* --- MAIN FDTD LOOP --- */
            my_ez[0] += CC*(my_hy[0] - my_hy[-MESH_SIZE] - my_hx[0] + my_hx[-1]);

            // pulse centered at PULSE_X, PULSE_Y with PULSE_WIDTH width and PULSE_SPREAD "height"
            // adds pulse every timestep at a constant intensity
            if (y == PULSE_Y && x >= MAX(PULSE_X-PULSE_WIDTH/2, 0) && x < MIN(PULSE_X+PULSE_WIDTH/2, MESH_SIZE)) {
                float pulse = (1.f / PULSE_SPREAD);
                pulse = exp(-0.5f * pulse * pulse) * (2.f / PULSE_WIDTH);
                my_ez[0] += pulse;
            }
        }

        __syncthreads();

        bool printout = --steps_until_printout == 0;
        for (unsigned int i = i0; i < (MESH_SIZE-2)*(MESH_SIZE-2); i += blockDim.x) {
            unsigned int x = i % MESH_SIZE, y = i / MESH_SIZE;
            if (x == 0 || y == 0 || x == MESH_SIZE - 1 || y == MESH_SIZE - 1) { continue; }
            float* my_ez = &ez[y*MESH_SIZE + x];
            float* my_hy = &hy[y*MESH_SIZE + x];
            float* my_hx = &hx[y*MESH_SIZE + x];

            // Calculate magnetic field in the X/Y
            my_hx[0] += CC*(my_ez[0] - my_ez[1]);
            my_hy[0] += CC*(my_ez[MESH_SIZE] - my_ez[0]);

            /* --- END OF MAIN FDTD LOOP --- */
            if (printout) {
                output[t/SAVE_EVERY_N_STEPS*MESH_SIZE_SQUARED + y*MESH_SIZE + x] = my_ez[0];
            }
        }
        if (printout) { steps_until_printout = SAVE_EVERY_N_STEPS; }

        __syncthreads();
    }
}

int main(int argc, const char *argv[]) {
    // parse args
    const char* output_loc = "./output-fdtd-2d-additive-cuda.npy";;
    if (argc > 1) {
        output_loc = argv[1];
    }
    const size_t num_elems = MESH_SIZE*MESH_SIZE, nbytes = num_elems * sizeof(float);
    float *d_ez, *d_hy, *d_hx, *d_output;
    // round up
    const size_t num_outputs = (NUM_TIMESTEPS+SAVE_EVERY_N_STEPS-1)/SAVE_EVERY_N_STEPS;
    Matrix *output = matrix_create_raw(num_outputs, num_elems);
    CHECK(cudaMalloc(&d_ez, nbytes));
    CHECK(cudaMalloc(&d_hy, nbytes));
    CHECK(cudaMalloc(&d_hx, nbytes));
    CHECK(cudaMalloc(&d_output, nbytes*num_outputs));
    CHECK(cudaMemset(d_ez, 0, nbytes));
    CHECK(cudaMemset(d_hy, 0, nbytes));
    CHECK(cudaMemset(d_hx, 0, nbytes));
    CHECK(cudaMemset(d_output, 0, nbytes*num_outputs));

    // start the clock
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    main_device<<<1, BLOCK_SIZE>>>(d_ez, d_hy, d_hx, d_output);
    CHECK(cudaDeviceSynchronize());

    // get the end and computation time
    clock_gettime(CLOCK_MONOTONIC, &end);
    double time = get_time_diff(&start, &end);
    printf("%f secs\n", time);

    CHECK(cudaMemcpy(output->data, d_output, num_outputs*nbytes, cudaMemcpyDeviceToHost));

    // save results
    matrix_to_npy_path(output_loc, output);
    
	//fclose(fp);
    matrix_free(output);
    CHECK(cudaFree(d_ez));
    CHECK(cudaFree(d_hy));
    CHECK(cudaFree(d_hx));
    CHECK(cudaFree(d_output));
    CHECK(cudaDeviceReset());
    return 0;
}