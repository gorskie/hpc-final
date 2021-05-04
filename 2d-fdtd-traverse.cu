/**
 * Proof of concept for the 2 dimensional FDTD simulation of EM waves
 * 
 * To compile the program:
 * gcc -Wall -O3 -march=native -c matrix.c util.c   # only need to run once
 * nvcc -arch=sm_20 -O3 *.cu *.o -o 2d-fdtd-transvers-cuda -lm
 * 
 * To run the program:
 *  ./2d-fdtd-transvers-cuda output.npy
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cuda_runtime.h>

#ifndef MESH_SIZE
#define MESH_SIZE 60
#endif
#ifndef NUM_TIMESTEPS
#define NUM_TIMESTEPS 100
#endif
#ifndef T0
#define T0 20 // >= 2
#endif
#ifndef PULSE_SPREAD
#define PULSE_SPREAD 6.f // gaussian
#endif
#ifndef PULSE_WIDTH
#define PULSE_WIDTH 40
#endif
#ifndef PULSE_Y
#define PULSE_Y 2
#endif
#ifndef PULSE_X
#define PULSE_X MESH_SIZE/2
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
void main_device(float* ez_row, float* hy_row, float* hx_row) {
    unsigned int ix = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int iy = blockIdx.y * blockDim.y + threadIdx.y;

    ez_row += ix*MESH_SIZE; hy_row += ix*MESH_SIZE; hx_row += ix*MESH_SIZE;

    hx_row[iy] += CC*(ez_row[iy] - ez_row[iy+1]);
    hy_row[iy] += CC*(ez_row[MESH_SIZE+iy] - ez_row[iy]);
}

int main() {
    size_t nbytes = MESH_SIZE*MESH_SIZE;
	float ez[nbytes], hx[nbytes], hy[nbytes];
    memset(ez, 0, nbytes); memset(hx, 0, nbytes); memset(hy, 0, nbytes);
    float *ez_row, *hy_row, *hx_row;
    float *d_ez_row, *d_hy_row, *d_hx_row;

    CHECK(cudaMalloc(&d_ez_row, nbytes));
    CHECK(cudaMalloc(&d_hy_row, nbytes));
    CHECK(cudaMalloc(&d_hx_row, nbytes));

    FILE* fp = fopen("2d_data.csv", "w");

    for (int t = 0; t < NUM_TIMESTEPS; t++) {

        /* --- MAIN FDTD LOOP --- */
        ez_row = ez; hy_row = hy; hx_row = hx;
        for (int i = 1; i < MESH_SIZE; i++) {
            ez_row += MESH_SIZE; hy_row += MESH_SIZE; hx_row += MESH_SIZE;
            for (int j = 1; j < MESH_SIZE; j++) {
                ez_row[j] += CC*(hy_row[j] - hy_row[j-MESH_SIZE] - hx_row[j] + hx_row[j-1]);
            }
        }

        // pulse centered at PULSE_X, PULSE_Y with PULSE_WIDTH width and PULSE_SPREAD "height"
        // peaks in intentsity at T0 and decays with time
        float pulse = (T0-t)*(1.f / PULSE_SPREAD);
        pulse = exp(-0.5f * pulse * pulse) * (2.f / PULSE_WIDTH);
        ez_row = ez+PULSE_Y*MESH_SIZE;
        for (int j = MAX(PULSE_X-PULSE_WIDTH/2, 0); j < MIN(PULSE_X+PULSE_WIDTH/2, MESH_SIZE); j++) {
            ez_row[j] = pulse;
        }

        /* --- boundary conditions --- */
        memset(ez, 0, MESH_SIZE*sizeof(float));
        memset(ez+MESH_SIZE*(MESH_SIZE-1), 0, MESH_SIZE*sizeof(float));
        for (int i = 1; i < MESH_SIZE-1; i++) { ez[i*MESH_SIZE] = ez[i*MESH_SIZE+MESH_SIZE-1] = 0.f; }

        CHECK(cudaMemcpy(d_ez_row, ez_row, nbytes, cudaMemcpyHostToDevice));
        CHECK(cudaMemcpy(d_hy_row, hy_row, nbytes, cudaMemcpyHostToDevice));
        CHECK(cudaMemcpy(d_hx_row, hx_row, nbytes, cudaMemcpyHostToDevice));
    
        const int n_block_size = 32;
        const int m_block_size = 32;
        dim3 grid((MESH_SIZE + n_block_size - 1) / n_block_size, (MESH_SIZE + m_block_size - 1) / m_block_size);
        dim3 block(n_block_size, m_block_size);
    
        main_device<<<grid, block>>>(d_ez_row, d_hy_row, d_hx_row);

        CHECK(cudaMemcpy(ez_row, d_ez_row, nbytes, cudaMemcpyDeviceToHost));
        CHECK(cudaMemcpy(hy_row, d_hy_row, nbytes, cudaMemcpyDeviceToHost));
        CHECK(cudaMemcpy(hx_row, d_hx_row, nbytes, cudaMemcpyDeviceToHost));
        CHECK(cudaDeviceSynchronize());

        /* --- MAIN FDTD LOOP --- */
        
        for (int i = 0; i < MESH_SIZE; i++) {
            fprintf(fp, "%g", ez[i*MESH_SIZE]);
            for (int j = 1; j < MESH_SIZE; j++) {
                fprintf(fp, ",%g", ez[i*MESH_SIZE+j]);
            }
            fprintf(fp, "\n");
        }
    }
    CHECK(cudaDeviceReset());
    CHECK(cudaFree(d_ez_row)); CHECK(cudaFree(d_hy_row)); CHECK(cudaFree(d_hx_row));
	fclose(fp);
}
