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
#define MESH_SIZE 100
#endif
/*
#ifndef SOURCE_POSITION_X
#define SOURCE_POSITION_X MESH_SIZE/2 // start in the middle
#endif
#ifndef SOURCE_POSITION_Y
#define SOURCE_POSITION_Y MESH_SIZE/2 // start in the middle
#endif
*/
#ifndef NUM_TIMESTEPS
#define NUM_TIMESTEPS 200
#endif
#ifndef TIMESTEP_SAVE
#define TIMESTEP_SAVE 5 // 20 CANNOT BE < 2
#endif
#ifndef CELL_SIZE
#define CELL_SIZE 3.e1
#endif
//#ifndef DZ
//#define DZ 0.01
//#endif
#ifndef PULSE_SPREAD
#define PULSE_SPREAD 36 // gaussian
#endif
//#define DT DZ/(2.*CELL_SIZE)
#define CC 0.5//CELL_SIZE*DT/DZ

void save_pos(Matrix* E_field, Matrix* output, size_t num_saves) {
    float* restrict E_row = E_field->data;
    float* restrict output_data = output->data;
    //size_t output_current_roww = num_saves - 1;
    // the current row will be the time step
    size_t pos = num_saves*MESH_SIZE*MESH_SIZE;
    //printf("%zu\n", output_curr_row+1);
    //printf("%zu\n", pos);

    // TODO: FIXME
    // populate the current row in the output matrix
    for (size_t i = 0; i < MESH_SIZE*MESH_SIZE; i++, pos++)
    {   
        // TODO: remove this conditional later
        if (!isfinite(E_row[i])) {
            printf("%f\n", E_row[i]);
        }
        if (E_row[i] != 0) {// && MESH_SIZE*MESH_SIZE%i==0) {
            printf("%zu %g\n", pos, E_row[i]);
        }
        /*if (i!=0 && (MESH_SIZE*MESH_SIZE)%i==0) {
            printf("%zu %g\n", pos, E_row[i]);
        }
        }*/
        output_data[pos] = E_row[i];
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
    Matrix* output_matrix = matrix_zeros(MESH_SIZE*(NUM_TIMESTEPS/TIMESTEP_SAVE), MESH_SIZE);
    
    // print out every TIMESTEP_SAVE steps
    size_t steps_until_printout = TIMESTEP_SAVE;

    // start timesteps
    for (size_t t = 1, num_saves = 0; t <= NUM_TIMESTEPS; t++, steps_until_printout--) {
        printf("%zu\n", steps_until_printout);

        /* --- MAIN FDTD LOOP ---
        D field
			for ( j=1; j < NY; j++) {
				for ( i=1; i < NX; i++) {
					dz[i][j] = dz[i][j] + 0.5*(hy[i][j] - hy[i-1][j] - hx[i][j] + hx[i][j-1]);
				}
			}
        */
        size_t index;
        // TODO: convert to one loop
        for (size_t i = 1; i < MESH_SIZE; i++) {
            for (size_t j = 1; j < MESH_SIZE; j++) {
                //index = (i-1)*MESH_SIZE+(j-1); 
                index = i*MESH_SIZE+j;
                // TODO: change order in which do math
                D_field->data[index] = D_field->data[index] + 0.5*
                    (H_field_y->data[index-MESH_SIZE] - H_field_y->data[index] + (H_field_x->data[index-1] - H_field_x->data[index]));
                if (!isfinite(D_field->data[index])) { printf("%f %f %f %f %f\n", D_field->data[index],
                H_field_y->data[index-MESH_SIZE], H_field_y->data[index], H_field_x->data[index], H_field_x->data[index-1]); 
                if (isnan(D_field->data[index])) {
                    return 0;
                }
                }
            }
        }
        /*
        for (size_t i = 1; i < MESH_SIZE; i++) {
            D_field->data[i] = D_field->data[i] + 0.5*(H_field_y->data[]);
        }
        */
        /*
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
        */
        
        /*
        Original by Will Langford
        pulse = exp(-0.5*(pow((t0-T)/spread,2.0)));
        #define NX 60
        #define NY 60
        jc = NX/2;
        jc = 60/2 = 30 - 20 = 10, 30 + 20 = 50
        range is jc-20, jc + 20
        jc = NX / 2; (world size/2)
			for ( j = jc-20; j < jc+20; j++) {
				// dz[1][j] = pulse;
				dz[2][j] = pulse/20.;
			}

        */
        float t_delta_amt = TIMESTEP_SAVE-t;
        float t_delta_amt_p = (t_delta_amt / PULSE_SPREAD);
        float t_delta_amt_p_squared = t_delta_amt_p * t_delta_amt_p;
        size_t boundary_offset = MESH_SIZE/6; // MESH_SIZE/6;
        float pulse_amount = exp(-0.5*t_delta_amt_p_squared);
        if (!isfinite(pulse_amount)) {
            printf("%f\n", pulse_amount);
            return 0;
        }
        /*for (size_t j = MESH_SIZE/6; j < 5*MESH_SIZE/6; j++) {
            D_field->data[MESH_SIZE+j] = pulse_amount/20;
        }
        */
        // Apply the source
        // this one works but isn't exactly what i expect
        // TODO: confirm correctness with jeff
        for (size_t bound_index = boundary_offset; bound_index < boundary_offset*5; bound_index++) {
            D_field->data[2*MESH_SIZE+bound_index] = pulse_amount;
        }
        //size_t matrix_center = MESH_SIZE/2;
        /*
        for (size_t index = boundary_offset; index < 3*boundary_offset; index++) {
            D_field->data[2*(MESH_SIZE-1)+index] = pulse_amount; // why cant this be the first row in the matrix?
        }
        */

        // Boundary conditions
        /*
        Original
        for ( j=1; j < NY; j++) {
            dz[0][j] = 0.0;
            dz[1][j] = 0.0;
            dz[NX][j] = 0.0;
            dz[NX-1][j] = 0.0;
        }
        for ( i=1; i < NY; i++) {
            dz[i][0] = 0.0;
            dz[i][1] = 0.0;
            dz[i][NY] = 0.0;
            dz[i][NY-1] = 0.0;
        }

        Ours
        (Something not working)
        for (size_t dz_j = 1; dz_j < MESH_SIZE; dz_j++) {
            //D_field->data[dz_j] = 0;
            //D_field->data[MESH_SIZE + dz_j] = 0;
            //D_field->data[(MESH_SIZE*(MESH_SIZE-1)) + dz_j] = 0;
            D_field->data[(MESH_SIZE*(MESH_SIZE)) + dz_j] = 0;
        }
        
        for (size_t dz_i = 1; dz_i < MESH_SIZE; dz_i++) {
            //D_field->data[MESH_SIZE] = 0;
            //D_field->data[MESH_SIZE + dz_i] = 0;
            //D_field->data[(MESH_SIZE-1)*dz_i] = 0;
            D_field->data[(MESH_SIZE)*dz_i] = 0;
        }
        */


        // Calculate electric field intensities
        // TODO: convert to 1 loop
        size_t current_row_ex = 0;
        size_t index_ex;// = current_row_ex;
        // TODO: figure out if index assignments in outer loop are unnecessary
        // TODO: possibly can combine this with the D field loop, but might slow down cache
        for (size_t ex_i = 1; ex_i < MESH_SIZE; ex_i++, current_row_ex+=MESH_SIZE, index_ex=current_row_ex) {
            for (size_t ex_j = 1; ex_j < MESH_SIZE; ex_j++, index_ex++) {
                // Not sure why this 1 needs to be here
                // TODO: determine if can remove without affecting output
                E_field->data[index_ex] = D_field->data[index_ex];
            }
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
        /*
        // Converted to one loop
        for (size_t hx_i = 0; hx_i < MESH_SIZE-1; hx_i++) {
            H_field_x->data[hx_i] = H_field_x->data[hx_i] + CC*(E_field->data[hx_i] - E_field->data[hx_i+1]);
        }
        */

        /*
        for ( j=0; j < NY-1; j++) {
				for ( i=0; i < NX-1; i++) {
					hx[i][j] = hx[i][j] + 0.5*(ez[i][j] - ez[i][j+1]);
 				}
			}

			for ( j=0; j < NY-1; j++) {
				for ( i=0; i < NX-1; i++) {
					hy[i][j] = hy[i][j] + 0.5*(ez[i+1][j] - ez[i][j]);
				}
			}
        */
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
        if (steps_until_printout == 1) {
            //printf("saving at step %zu\n", t);
            save_pos(E_field, output_matrix, num_saves++);
            // reset the timestep print counter
            // note that this will still decrement at the end of the loop
            steps_until_printout += TIMESTEP_SAVE;
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

