#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

//#include <cooperative_groups.h>

#include "matrix.h"
#include "util.h"

#ifndef MESH_SIZE
#define MESH_SIZE 100
#endif
#define MESH_SIZE_SQUARED MESH_SIZE*MESH_SIZE
#ifndef NUM_TIMESTEPS
#define NUM_TIMESTEPS 100
#endif
#ifndef SAVE_EVERY_N_STEPS
#define SAVE_EVERY_N_STEPS 5
#endif
#ifndef T0
#define T0 20 // >= 2 This is the peak of when the source is applied
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
#ifndef N_BLOCK_SIZE
#define N_BLOCK_SIZE 32
#endif
#ifndef M_BLOCK_SIZE
#define M_BLOCK_SIZE 32
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
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;

    //cooperative_groups::grid_group grid = cooperative_groups::this_grid();

    float* my_ez = &ez[i*MESH_SIZE + j];
    float* my_hy = &hy[i*MESH_SIZE + j];
    float* my_hx = &hx[i*MESH_SIZE + j];
    my_ez[0] = 0;
    my_hy[0] = 0;
    my_hx[0] = 0;

    if (1 <= i && i < MESH_SIZE && 1 <= j && j < MESH_SIZE) {
        for (int t = 0, steps_until_printout = SAVE_EVERY_N_STEPS; t < NUM_TIMESTEPS; t++) {
            
            __syncthreads();

            /* --- MAIN FDTD LOOP --- */
            my_ez[0] += CC*(my_hy[0] - my_hy[-MESH_SIZE] - my_hx[0] + my_hx[-1]);

            // pulse centered at PULSE_X, PULSE_Y with PULSE_WIDTH width and PULSE_SPREAD "height"
            // peaks in intensity at T0 and decays with time
            if (i = PULSE_Y && j >= MAX(PULSE_X-PULSE_WIDTH/2, 0) && j < MIN(PULSE_X+PULSE_WIDTH/2, MESH_SIZE)) {
                float pulse = (T0-t)*(1.f / PULSE_SPREAD);
                pulse = exp(-0.5f * pulse * pulse) * (2.f / PULSE_WIDTH);
                my_ez[0] += pulse;
            }

            __syncthreads();

            // Calculate magnetic field in the X
            my_hx[0] += CC*(my_ez[0] - my_ez[1]);

            // Calculate magnetic field in the Y
            my_hy[0] += CC*(my_ez[MESH_SIZE] - my_ez[0]);

            /* --- END OF MAIN FDTD LOOP --- */
            if (--steps_until_printout == 0) {
                output[t/SAVE_EVERY_N_STEPS*MESH_SIZE_SQUARED + i*MESH_SIZE + j] = my_ez[0];
                steps_until_printout = SAVE_EVERY_N_STEPS;
            }
        }
    }
}

int main(int argc, const char *argv[]) {
    // parse args
    const char* output_loc = "./output-fdtd-2d-transverse.npy";
    if (argc > 1) {
        output_loc = argv[1];
    }

    const size_t num_elems = MESH_SIZE*MESH_SIZE, nbytes = num_elems * sizeof(float);
	float ez[num_elems], hx[num_elems], hy[num_elems];
    memset(ez, 0, nbytes); memset(hx, 0, nbytes); memset(hy, 0, nbytes);
    float *d_ez, *d_hy, *d_hx, *d_output;
    const size_t num_outputs = (NUM_TIMESTEPS+SAVE_EVERY_N_STEPS-1)/SAVE_EVERY_N_STEPS;
    Matrix *output = matrix_create_raw(num_outputs, num_elems);
    cudaStream_t cuda_stream;
    CHECK(cudaStreamCreate(&cuda_stream));

    CHECK(cudaMalloc(&d_ez, nbytes));
    CHECK(cudaMalloc(&d_hy, nbytes));
    CHECK(cudaMalloc(&d_hx, nbytes));
    CHECK(cudaMalloc(&d_output, nbytes*num_outputs));

    dim3 grid((MESH_SIZE + N_BLOCK_SIZE - 1) / N_BLOCK_SIZE,
              (MESH_SIZE + M_BLOCK_SIZE - 1) / M_BLOCK_SIZE);
    dim3 block(N_BLOCK_SIZE, M_BLOCK_SIZE);
    main_device<<<grid, block>>>(d_ez, d_hy, d_hx, d_output);
    void* args[] = {&d_ez, &d_hx, &d_hy, &d_output};
    CHECK(cudaDeviceSynchronize());

    CHECK(cudaMemcpy(output->data, d_output, num_outputs, cudaMemcpyDeviceToHost));
    
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
