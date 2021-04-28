/**
 * Proof of concept for the 2 dimensional FDTD simulation of EM waves
 * 
 * To compile the program:
 * gcc -Wall -g -O3 -march=native 2d-fdtd-s.c matrix.c util.c -o 2d-fdtd-s -lm
 * 
 * To compile the program with different global parameters:
 * gcc -Wall -g -O3 -march=native 2d-fdtd-s.c matrix.c util.c -o 2d-fdtd-s -lm -DMESH_SIZE 60 -DSOURCE_POSITION_X 30 -DSOURCE_POSITION_Y 30 -DNUM_TIMESTEPS 100 -DTIMESTEP_SAVE 20 -DCELL_SIZE 3.e8 -DZ 0.01 -DPULSE_SPREAD 10
s * 
 * To run the program:
 *  ./2d-fdtd-s output-2d-s.npy
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <time.h>

#include "matrix.h"
#include "util.h"


#ifndef MESH_SIZE // Not supporting non-square FDTD spaces
#define MESH_SIZE 50
#endif
#ifndef SOURCE_POSITION_X
#define SOURCE_POSITION_X MESH_SIZE/2 // start in the middle
#endif
#ifndef SOURCE_POSITION_Y
#define SOURCE_POSITION_Y MESH_SIZE/2 // start in the middle
#endif
#ifndef NUM_TIMESTEPS
#define NUM_TIMESTEPS 100
#endif
#ifndef TIMESTEP_SAVE
#define TIMESTEP_SAVE 20 // 20
#endif
#ifndef CELL_SIZE
#define CELL_SIZE 3.e8
#endif
#ifndef DZ
#define DZ 0.01
#endif
#ifndef PULSE_SPREAD
#define PULSE_SPREAD 8 // gaussian
#endif
#define DT DZ/(2.*CELL_SIZE)
#define CC CELL_SIZE*DT/DZ

// Structs
// TODO: determine if needed
/*
typedef struct _electric_field {
    float x[MESH_SIZE_X];
    float y[MESH_SIZE_Y];
} electric_field;

typedef struct _curl_field {
    float x[MESH_SIZE_X];
    float y[MESH_SIZE_Y];
} curl_field;
*/

// FIXME: output correctly to 2d
void save_pos(Matrix* E_field, Matrix* output, size_t output_curr_row) {
    float* restrict E_row = E_field->data;
    float* restrict output_data = output->data;

    // the current row will be the time step
    size_t pos = output_curr_row*MESH_SIZE*MESH_SIZE;
    //printf("%zu\n", output_curr_row+1);
    //printf("%zu\n", pos);

    // TODO: FIXME
    // populate the current row in the output matrix
    for (size_t i = 0; i < MESH_SIZE*MESH_SIZE; i++)
    {   
        // TODO: remove this conditional later
        if (E_row[i] != 0 && MESH_SIZE*MESH_SIZE%i==0) {
            printf("%zu %g\n", pos, E_row[i]);
        } 
        output_data[pos] = E_row[i];
        pos++;
    }
}


// TODO: finish converting to 2d
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

    // initializes the matrixes for the electric field intensities at x, y and the magnetic fields at x, y components
    Matrix* D_field = matrix_zeros(MESH_SIZE, MESH_SIZE);
    Matrix* E_field = matrix_zeros(MESH_SIZE, MESH_SIZE);
    Matrix* H_field_x = matrix_zeros(MESH_SIZE, MESH_SIZE);
    Matrix* H_field_y = matrix_zeros(MESH_SIZE, MESH_SIZE);
    Matrix* output_matrix = matrix_zeros(MESH_SIZE*NUM_TIMESTEPS, MESH_SIZE);
    
    // print out every TIMESTEP_SAVE steps
    size_t steps_until_printout = TIMESTEP_SAVE-1;

    // start timesteps
    for (size_t t = 0; t < NUM_TIMESTEPS; t++) {
        //printf("t = %zu\n", t);
        //printf("steps until printout: %zu\n", steps_until_printout);
        // calculate the intensities of the electric field at nodes
        // is one step ahead of the curl (h field)
        
        // calculate D field
        size_t prev_row_dz = 0;
        size_t current_row_dz = MESH_SIZE;
        size_t prev_row_index_dz = prev_row_dz;
        size_t index_dz = current_row_dz;
        // TODO: figure out if index assignments in outer loop are unnecessary
        for (size_t dz_i = 1; dz_i < MESH_SIZE; dz_i++, prev_row_dz+=MESH_SIZE, current_row_dz+=MESH_SIZE, prev_row_index_dz=prev_row_dz, index_dz=current_row_dz) {
            for (size_t dz_j = 1; dz_j < MESH_SIZE; dz_j++, prev_row_index_dz++, index_dz++) {
                D_field->data[index_dz] = D_field->data[index_dz] + CC * (H_field_y->data[index_dz] - H_field_y->data[prev_row_index_dz] -
                                                        H_field_x->data[index_dz] + H_field_y->data[index_dz-1]);
            }
            //prev_row_dz += MESH_SIZE;
            //current_row_dz += MESH_SIZE;
        }
        
        // apply source
        float t_delta_amt = (t-(TIMESTEP_SAVE-1));
        float t_delta_amt_p = (t_delta_amt / PULSE_SPREAD);
        float t_delta_amt_p_squared = t_delta_amt_p * t_delta_amt_p;
        //printf("dt = %f\np = %d\n", t_delta_amt, PULSE_SPREAD);
        //printf("e^(-1/2 * (dt/p)^2) = %f\n", exp(-0.5*t_delta_amt_p_squared));
        // TODO: reduce math, move this elsewhere
        size_t center = MESH_SIZE*MESH_SIZE/2 + MESH_SIZE/2;
        D_field->data[center] = exp(-0.5*t_delta_amt_p_squared);

        // Calculate electric field intensities
        // TODO: convert to 1 loop
        size_t current_row_ex = 0;
        size_t index_ex;// = current_row_ex;
        // TODO: figure out if index assignments in outer loop are unnecessary
        // TODO: possibly can combine this with the D field loop, but might slow down cache
        for (size_t ex_i = 1; ex_i < MESH_SIZE; ex_i++, current_row_ex+=MESH_SIZE, index_ex=current_row_ex) {
            for (size_t ex_j = 1; ex_j < MESH_SIZE; ex_j++, index_ex++) {
                //size_t index = current_row_ex+ex_j;
                //H_field_x->data[] = H_field->data[k] + CC_X*(E_field->data[k]-H_field->data[k+1]);
                // Not sure why this 1 needs to be here
                // TODO: determine if can remove without affecting output
                E_field->data[index_ex] = 1*D_field->data[index_ex];
            }
            //current_row_ex += MESH_SIZE;
        }

        // Calculate magnetic field in the X
        //size_t prev_row_hx = 0;
        // TODO: figure out if index assignments in outer loop are unnecessary
        size_t current_row_hx = MESH_SIZE;
        size_t index_hx = 0;
        for (size_t hx_i = 0; hx_i < MESH_SIZE-1; hx_i++, current_row_hx+=MESH_SIZE, index_hx=current_row_hx) {
            for (size_t hx_j = 0; hx_j < MESH_SIZE-1; hx_j++, index_hx++) {
                H_field_x->data[index_hx] = H_field_x->data[index_hx] + CC*(E_field->data[index_hx] - E_field->data[index_hx+1]);
            }
        }

        // Calculate magnetic field in the Y
        size_t current_row_hy = 0;
        size_t next_row_hy = MESH_SIZE;
        size_t index_hy = current_row_hy;
        size_t next_row_index_hy = next_row_hy;
        // TODO: figure out if index assignments in outer loop are unnecessary
        for (size_t hy_i = 0; hy_i < MESH_SIZE-1; hy_i++, current_row_hy+=MESH_SIZE, next_row_hy+=MESH_SIZE, index_hy=current_row_hy, next_row_index_hy=next_row_hy) {
            for (size_t hy_j = 0; hy_j < MESH_SIZE-1; hy_j++, index_hy++, next_row_index_hy++) {
                //hy[i][j] = hy[i][j] + 0.5*(ez[i+1][j] - ez[i][j]);
                H_field_y->data[index_hy] = H_field_y->data[index_hy] + CC*(E_field->data[next_row_index_hy] - E_field->data[index_hy]);
            }
        }

        // save values at the timestep
        if (steps_until_printout-1 == 0) {
            printf("saving at step %zu\n", t);
        }
        if (steps_until_printout-- == 0) {
            //printf("actual current timestep before saving: %zu\n", t);
            //printf("value being sent to save: (t+1)/TIMESTEP_SAVE = %zu/%d = %zu\n", t+1, TIMESTEP_SAVE, (t+1)/TIMESTEP_SAVE);
            // the current row is the human-readable timestep divided by the frequency of timestep saves
            save_pos(E_field, output_matrix, (t+1)/TIMESTEP_SAVE);
            // reset the timestep print counter
            steps_until_printout = TIMESTEP_SAVE - 1;
        }
    }

    // save results
    matrix_to_npy_path(output_loc, output_matrix);
    
    // free the matrices
    matrix_free(D_field);
    matrix_free(E_field);
    matrix_free(H_field_x);
    matrix_free(H_field_y);
    matrix_free(output_matrix);

    return 0;
}

