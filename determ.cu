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
#include "gputimer.hpp"

Options options;
const Options& opts = options;

float* init_f(const VelocityGrid& h_vgrid)
{
	float * f = new float[opts.nxyz*opts.npx];
	float u1 = 1.75;
	float v1 = 0.25;

	for (int i = 0, index = 0; i < opts.nx; i++) 
		for (int j = 0; j < opts.ny; j++) 
			for (int k = 0; k < opts.nz; k++, index++) {
				float u = h_vgrid.u[i];
				float v = h_vgrid.v[j];
				float w = h_vgrid.w[k];
	
				for (int ix = 0; ix < opts.npx; ix++) 
					f[index * opts.npx + ix] = exp (-(u-u1)*(u-u1)  - (v-v1)*(v-v1) - w*w); 
			}

	return f;
}

int main()
{
	load_options(options, "gpu.conf");
	print_options(options);

	/* VELOCITY GRID INITIALIZATION */
	const Vector3i n_points(opts.nx, opts.ny, opts.nz);
	const Vector3f v_min(-5.0, -5.0, -5.00);
	const Vector3f v_max(5.0, 5.0, 5.0);
	const float R = 16.5; // RADIUS OF SPHERE
	/*  VELOCITY GRID INITIALIZATION */

	VelocityGrid h_vgrid(n_points, v_min, v_max, R);
	VelocityGrid d_vgrid = h_vgrid.device_clone();

	float * correction_array;
	init_correction_array (&correction_array, d_vgrid, opts);

	float * a;
	float * b;
  	init_matrices (d_vgrid, &b, &a, opts);

  float * h_f = init_f(h_vgrid);//new float[opts.nxyz*opts.npx];
  float * h_inverse_integral = new float[opts.nxyz*opts.npx];
  float * h_direct_integral = new float[opts.nxyz*opts.npx];
  float * d_f;
  float * d_inverse_integral;
  float * d_direct_integral;
  int index, step, out_step = 5;

  cudaEvent_t start, stop;
  float elapsedTime;

  int i, j, k;
  float h_time = 0;
  float * d_time;


  int ix;


  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaMalloc ((void **) &d_f, opts.nxyz * opts.npx * sizeof (float));
  cudaMalloc ((void **) &d_direct_integral, opts.nxyz * opts.npx * sizeof (float));
  cudaMalloc ((void **) &d_inverse_integral, opts.nxyz * opts.npx * sizeof (float));
  cudaMalloc ((void **) &d_time,  sizeof (float));



  /* INITIALIZATION OF THE DISTRIBUTION FUNCTION h_f */ 
  /*              h_f = F(u,v,w)                     */ 
  /*
	for (i = 0, index = 0; i < opts.nx; i++) 
		for (j = 0; j < opts.ny; j++) 
			for (k = 0; k < opts.nz; k++, index++) {
				float u = h_vgrid.u[i];
				float v = h_vgrid.v[j];
				float w = h_vgrid.w[k];
	
				for (ix = 0; ix < opts.npx; ix++) 
					h_f[index * opts.npx + ix] = exp (-(u-u1)*(u-u1)  - (v-v1)*(v-v1) - w*w); 
			}
			*/
  
  /* INITIALIZATION OF THE DISTRIBUTION FUNCTION h_f */ 
  
  

	cudaMemcpy (d_f, h_f, opts.nxyz * opts.npx * sizeof (float), cudaMemcpyHostToDevice);
	cudaMemcpy (d_time, &h_time, sizeof (float), cudaMemcpyHostToDevice);
  
  
	cudaThreadSynchronize ();
	cudaEventRecord(start, 0);
#if 1  
  
	for (step = 0; step <= 0; step ++) 
	{
		cudaMemcpy (&h_time, d_time, sizeof (float), cudaMemcpyDeviceToHost);
		printf ("STEP=%d, TIME=%f\n", step, time);
		if (step % out_step == 0) 
			cudaMemcpy (h_f, d_f, opts.nxyz * opts.npx * sizeof (float), cudaMemcpyDeviceToHost);

		evolve_colisions (d_f, d_direct_integral, d_inverse_integral, b, a, correction_array, opts);
		relax_f (d_f, d_direct_integral, d_inverse_integral, 0.5, d_time, opts);
	}
#endif
  
	cudaThreadSynchronize ();

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);

	printf ("COLLISIONS EVOLVED: %s\n", cudaGetErrorString (cudaGetLastError ()));

	cudaMemcpy (h_direct_integral, d_direct_integral, opts.nxyz * sizeof (float), cudaMemcpyDeviceToHost);
	cudaMemcpy (h_inverse_integral, d_inverse_integral, opts.nxyz * sizeof (float), cudaMemcpyDeviceToHost);

	std::ofstream fresults("results.out");

#if 1
  /* PRINTING THE RESULTS OF COMPUTATIONS */ 
  for (index = 0; index < opts.nxyz; index ++) 
    //    if (h_f[index]>1.e-1)
    if (h_vgrid.u_index[index] < h_vgrid.n_u/2 && 
	h_vgrid.v_index[index] < h_vgrid.n_v/2 &&
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
	float v = h_vgrid.v[j];//V_1 (h_vgrid.float_params, h_vgrid.int_params)[j];
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

  free_device_velocity_grid (d_vgrid);
  free_host_velocity_grid (h_vgrid);

	printf ("END_PROGRAM: %s\n", 
	cudaGetErrorString (cudaGetLastError ()));
}
