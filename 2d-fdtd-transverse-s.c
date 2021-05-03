#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

int main() {
	float ez[MESH_SIZE*MESH_SIZE] = {0}, hx[MESH_SIZE*MESH_SIZE] = {0}, hy[MESH_SIZE*MESH_SIZE] = {0};
    float *ez_row, *hy_row, *hx_row;

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

        ez_row = ez; hx_row = hx;
        for (int i = 0; i < MESH_SIZE-1; i++) {
            for (int j = 0; j < MESH_SIZE-1; j++) {
                hx_row[j] += CC*(ez_row[j] - ez_row[j+1]);
            }
            ez_row += MESH_SIZE; hx_row += MESH_SIZE;
        }

        ez_row = ez; hy_row = hy;
        for (int i = 0; i < MESH_SIZE-1; i++) {
            for (int j = 0; j < MESH_SIZE-1; j++) {
                hy_row[j] += CC*(ez_row[MESH_SIZE+j] - ez_row[j]);
            }
            ez_row += MESH_SIZE; hy_row += MESH_SIZE;
        }

        /* --- MAIN FDTD LOOP --- */
        
        for (int i = 0; i < MESH_SIZE; i++) {
            fprintf(fp, "%g", ez[i*MESH_SIZE]);
            for (int j = 1; j < MESH_SIZE; j++) {
                fprintf(fp, ",%g", ez[i*MESH_SIZE+j]);
            }
            fprintf(fp, "\n");
        }
    }
	fclose(fp);
}
