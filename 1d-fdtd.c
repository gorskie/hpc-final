/**
 * Proof of concept for the 1 dimensional FDTD simulation of EM waves
 * 
 * TODO: fix compile statement
 * To compile the program:
 * gcc -Wall -O3 -march=native 1d-fdtd.c matrix.c util.c -o 1d-fdtd -lpthread -lm
 * 
 * To compile the program with different global parameters:
 * gcc -Wall -O3 -march=native 1d-fdtd.c matrix.c util.c -o 1d-fdtd -lpthread -lm -DMESH_SIZE 60 -DSOURCE_POSITION 30 -DNUM_TIMESTEPS 100 -DTIMESTEP_SAVE 20 -DCELL_SIZE 3.e8 -DX 0.01 -DPULSE_SPREAD 10
 * TODO: fix this
 * 
 * To run the program:
 *  ./1d-fdtd output.npy
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <time.h>

#include "matrix.h"
#include "util.h"


#ifndef MESH_SIZE
#define MESH_SIZE 50 // size of fdtd space
#endif
#ifndef SOURCE_POSITION
#define SOURCE_POSITION MESH_SIZE/2 // start in the middle
#endif
#ifndef NUM_TIMESTEPS
#define NUM_TIMESTEPS 100
#endif
#ifndef TIMESTEP_SAVE
#define TIMESTEP_SAVE 20
#endif
#ifndef CELL_SIZE
#define CELL_SIZE 3.e8
#endif
#ifndef DX
#define DX 0.01
#endif
#ifndef PULSE_SPREAD
#define PULSE_SPREAD 8 // gaussian
#endif

#define DT DX/(2.*CELL_SIZE)
#define CC CELL_SIZE*DT/DX

void save_pos(Matrix* E_field, Matrix* output, size_t output_curr_row) {
    float* restrict E_row = &E_field->data[output_curr_row];
    float* restrict output_data = (output)->data;

    // the current row will be the time step
    size_t pos = output_curr_row*MESH_SIZE;
    // populate the current row in the output matrix
    for (size_t i = 0; i < MESH_SIZE; i++)
    {
        output_data[pos] = E_row[pos]; 
        pos++;
    }
}

int main(int argc, const char* argv[]) {
    // parse args
    const char* output_loc;
    if (argc > 1) {
        output_loc = argv[1];
    }
    else {
        output_loc = "./output.npy";
    }
    printf("%s\n", output_loc);
    // initializes the matrixes for the electric field and curl field
    Matrix* E_field = matrix_zeros(1, MESH_SIZE);
    Matrix* H_field = matrix_zeros(1, MESH_SIZE);
    Matrix* output_matrix = matrix_zeros(NUM_TIMESTEPS, MESH_SIZE);
    
    // print out every TIMESTEP_SAVE steps
    size_t steps_until_printout = TIMESTEP_SAVE;

    // start timesteps
    for (size_t t = 0; t < NUM_TIMESTEPS; t++) {
        // calculate the intensities of the electric field at nodes
        // is one step ahead of the curl (h field)
        
        // TODO: check boundaries
        for (size_t k = 1; k < MESH_SIZE; k++) {
            // potentially store the delta H_field from k-1 to k outside of this loop at each step?
            E_field->data[k] = E_field->data[k] + CC * (H_field->data[k-1]-H_field->data[k]);
        }
        
        // apply source
        float t_delta_amt = (t-TIMESTEP_SAVE) * (t-TIMESTEP_SAVE);
        E_field->data[SOURCE_POSITION] = exp(-0.5*(t_delta_amt/PULSE_SPREAD));

        // curl field
        // one step behind the electric field
        for (size_t k = 0; k < MESH_SIZE; k++) {
            H_field->data[k] = H_field->data[k] + CC*(E_field->data[k]-H_field->data[k+1]);
        }

        // save values at the timestep
        if (steps_until_printout-- == 0) {
            save_pos(E_field, output_matrix, t);
            // reset the timestep print counter
            steps_until_printout = TIMESTEP_SAVE;
        }
    }

    // save results
    matrix_to_npy_path(output_loc, output_matrix);
    
    // free the matrices
    matrix_free(E_field);
    matrix_free(H_field);
    matrix_free(output_matrix);
    return 0;
}

