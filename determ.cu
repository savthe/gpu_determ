/*
 * determ.cu
 *
 *  Created on: 17.02.2010
 *      Author: serge
 */

#include <fstream>
#include <cuda.h>
#include <stdio.h>
#include <assert.h>

#include "velocity.h"
#include "compute_matrix.h"
#include "collision_integral.h"
#include "options.hpp"
#include "velocity_grid.hpp"

Options options;
const Options& opts = options;

int main()
{
	load_options(options, "gpu.conf");
	print_options(options);
	VelocityGrid h_vgrid, d_vgrid;

  /* VELOCITY GRID INITIALIZATION */
  IntVector n_points = {opts.nx, opts.ny, opts.nz};
  DVector v_min = {-5.0, -5.0, -5.0};
  DVector v_max = {5.0, 5.0, 5.0};
  float R = 16.5; // RADIUS OF SPHERE
  /*  VELOCITY GRID INITIALIZATION */

  float * b;
  float * a;
  float * h_f;
  float * d_f;
  float * h_inverse_integral;
  float * d_inverse_integral;
  float * h_direct_integral;
  float * d_direct_integral;
  float * correction_array;
  int index, step, out_step = 5;

  cudaEvent_t start, stop;
  float elapsedTime;

  int i, j, k;
  float h_time = 0;
  float * d_time;

  float u1 = 1.75;
  float v1 = 0.25;

  int ix;


  //init_velocity_grid (h_vgrid, n_points, v_min, v_max, R);
  h_vgrid.init(n_points, v_min, v_max, R);


  init_device_velocity_grid (h_vgrid, d_vgrid);

  init_correction_array (&correction_array, d_vgrid, opts);

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);

  init_matrices (d_vgrid, &b, &a, opts);
  
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);

  printf ("MATRICES INITIALIZATION TOOK %f MS!\n", elapsedTime);


  cudaMalloc ((void **) &d_f, 
	      opts.nxyz * opts.npx * sizeof (float));
  cudaMalloc ((void **) &d_direct_integral, 
	      opts.nxyz * opts.npx * sizeof (float));
  cudaMalloc ((void **) &d_inverse_integral, 
	      opts.nxyz * opts.npx * sizeof (float));
  cudaMalloc ((void **) &d_time,  sizeof (float));

  h_f = (float *) malloc (opts.nxyz * opts.npx * sizeof (float));
  h_direct_integral = (float *) malloc (opts.nxyz * opts.npx * sizeof (float));
  h_inverse_integral = (float *) malloc (opts.nxyz * opts.npx * sizeof (float));


  /* INITIALIZATION OF THE DISTRIBUTION FUNCTION h_f */ 
  /*              h_f = F(u,v,w)                     */ 
  for (i = 0, index = 0; i < N_X; i++) 
    for (j = 0; j < N_Y; j++) 
      for (k = 0; k < N_Z; k++, index++) {
	
	float u = h_vgrid.u[i];//U_1 (h_vgrid.float_params, h_vgrid.int_params)[i];
	float v = h_vgrid.v[j];//V_1 (h_vgrid.float_params, h_vgrid.int_params)[j];
	float w = h_vgrid.w[k];//W_1 (h_vgrid.float_params, h_vgrid.int_params)[k];
	
	for (ix = 0; ix < opts.npx; ix++) {
	  h_f[index * opts.npx + ix] = exp (-(u-u1)*(u-u1)  - (v-v1)*(v-v1) - w*w); 
	  //	  h_f[index * opts.npx + ix] += exp (-(u-u2)*(u-u2)  - (v-v2)*(v-v2) - w*w); 
	}
      }
  
  /* INITIALIZATION OF THE DISTRIBUTION FUNCTION h_f */ 
  
  

  cudaMemcpy (d_f, h_f, opts.nxyz * opts.npx * sizeof (float),
	      cudaMemcpyHostToDevice);
  cudaMemcpy (d_time, &h_time, sizeof (float),
	      cudaMemcpyHostToDevice);
  
  
  cudaThreadSynchronize ();
  cudaEventRecord(start, 0);
#if 1  
  
  for (step = 0; step <= 0; step ++) {
    cudaMemcpy (&h_time, d_time, sizeof (float),
		cudaMemcpyDeviceToHost);
    printf ("STEP=%d, TIME=%f\n", step, time);
    if (step % out_step == 0) {
      cudaMemcpy (h_f, d_f, 
		  opts.nxyz * opts.npx * sizeof (float),
		  cudaMemcpyDeviceToHost);
      /*
      sprintf (filename_x, "RES1/x_data_%03d.dat", step);
      sprintf (filename_y, "RES1/y_data_%03d.dat", step);
      sprintf (filename_z, "RES1/z_data_%03d.dat", step);

      fp_x = fopen (filename_x, "w+");
      fprintf (fp_x, "T=%f\n", h_time);
      for (i = 0; i < N_X; i++) {
	fprintf (fp_x, "%g", U_1 (h_float_params, h_int_params)[i]);
	for (j = 0; j < N_Y; j++)
	  for (k = 0; k < N_Z; k++) {
	    fprintf (fp_x, " %g", h_f[i * N_YZ + j * N_Z + k]);
	  }
	fprintf (fp_x, "\n");
      }
      fclose (fp_x);

      fp_y = fopen (filename_y, "w+");
      fprintf (fp_y, "T=%f\n", h_time);
      for (j = 0; j < N_Y; j++) {
	fprintf (fp_y, "%g", V_1 (h_float_params, h_int_params)[j]);
	for (i = 0; i < N_X; i++)
	  for (k = 0; k < N_Z; k++) {
	    fprintf (fp_y, " %g", h_f[i * N_YZ + j * N_Z + k]);
	  }
	fprintf (fp_y, "\n");
      }
      fclose (fp_y);

      fp_z = fopen (filename_z, "w+");
      fprintf (fp_z, "T=%f\n", h_time);
      for (k = 0; k < N_Z; k++) {
	fprintf (fp_z, "%g", W_1 (h_float_params, h_int_params)[k]);
	for (j = 0; j < N_Y; j++)
	  for (i = 0; i < N_X; i++) {
	    fprintf (fp_z, " %g", h_f[i * N_YZ + j * N_Z + k]);
	  }
	fprintf (fp_z, "\n");
      }
      fclose (fp_z);
      */
    }
    evolve_colisions (d_f, d_direct_integral, d_inverse_integral,
		      b, a, correction_array, opts);
    relax_f (d_f, d_direct_integral, d_inverse_integral, 0.5,
	     d_time, opts);
  }
#endif
  
  cudaThreadSynchronize ();

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  printf ("COLLISIONS EVOLVED: %s\n", 
	  cudaGetErrorString (cudaGetLastError ()));

  cudaMemcpy (h_direct_integral, d_direct_integral, 
	      opts.nxyz * sizeof (float),
	      cudaMemcpyDeviceToHost);

  cudaMemcpy (h_inverse_integral, d_inverse_integral, 
	      opts.nxyz * sizeof (float),
	      cudaMemcpyDeviceToHost);

  std::ofstream fresults("results.out");

#if 1
  /* PRINTING THE RESULTS OF COMPUTATIONS */ 
  for (index = 0; index < opts.nxyz; index ++) 
    //    if (h_f[index]>1.e-1)
    if (U_INDEX (h_vgrid.int_params)[index] < h_vgrid.n_u/2 && 
	V_INDEX (h_vgrid.int_params)[index] < h_vgrid.n_v/2 &&
	h_vgrid.w_index[index] < h_vgrid.n_w/2)
	//W_INDEX (h_vgrid.int_params)[index] < h_vgrid.n_w/2)
    //  printf ("INDEX=%d, I=%d, J=%d, K=%d, F=%e, DC=%e, IC=%e, CI=%e, ACCURACY=%e\n", 
	fresults << index << ' ' // index
		<< index / N_YZ << ' ' // I
		<< (index % N_YZ) / N_Z << ' ' // J
	    << (index % N_YZ) % N_Z << ' ' // K
		<< h_f[index] << ' ' // F
		<< h_direct_integral[index] << ' ' // DC
		<< h_inverse_integral[index] << ' ' // IC
	    << - h_f[index]*h_direct_integral[index] + h_inverse_integral[index] << ' ' // CI
	    << -h_f[index]*h_direct_integral[index]/h_inverse_integral[index]+1.0 << std::endl; //ACCURACY
  /* PRINTING THE RESULTS OF COMPUTATIONS */ 
#endif
  fresults.close();

  /* EVALUATION OF MOMENTS OF COILLISION INTEGRAL */ 
  float mom0 = 0, momU = 0, momV = 0, momW = 0, mom2 = 0;
  for (i = 0, index = 0; i < N_X; i++) 
    for (j = 0; j < N_Y; j++)
      for (k = 0; k < N_Z; k++, index++) {
	float ci = - h_f[index] * h_direct_integral[index] + 
	  h_inverse_integral[index];
	float u = h_vgrid.u[i];//U_1 (h_vgrid.float_params, h_vgrid.int_params)[i];
	//float v = V_1 (h_vgrid.float_params, h_vgrid.int_params)[j];
	float v = h_vgrid.v[j];//V_1 (h_vgrid.float_params, h_vgrid.int_params)[j];
	//float w = W_1 (h_vgrid.float_params, h_vgrid.int_params)[k];
	float w = h_vgrid.w[k];
	mom0 += ci * h_vgrid.d3v;
	momU += ci * u * h_vgrid.d3v;
	momV += ci * v * h_vgrid.d3v;
	momW += ci * w * h_vgrid.d3v;
	mom2 += ci * (u*u + v*v + w*w) * h_vgrid.d3v;
      }

  printf ("M0=%e, MU=%e, MV=%e, MW=%e, M2=%e\n", mom0, momU,
	  momV, momW, mom2);
  /* EVALUATION OF MOMENTS OF COILLISION INTEGRAL */ 

  cudaEventElapsedTime(&elapsedTime, start, stop);

  printf ("COLLISION INTEGRAL EVALUATION TOOK %f MS!\n", elapsedTime);

  cudaEventDestroy (start);
  cudaEventDestroy (stop);

  cudaFree (d_f);
  cudaFree (d_direct_integral);
  cudaFree (d_inverse_integral);
  cudaFree (correction_array);

  free (h_f);
  free (h_direct_integral);
  free (h_inverse_integral);

  cudaFree (b);
  cudaFree (a);

  /*
  int * h_int_params = h_vgrid.int_params;
  float * h_float_params = h_vgrid.float_params;
  int * d_int_params = d_vgrid.int_params;
  float * d_float_params = d_vgrid.float_params;
  */
  free_device_velocity_grid (d_vgrid);
  free_host_velocity_grid (h_vgrid);
//  free_device_velocity_grid (&d_int_params, &d_float_params);

 // free_host_velocity_grid (&h_int_params, &h_float_params);

  printf ("END_PROGRAM: %s\n", 
	  cudaGetErrorString (cudaGetLastError ()));
}
