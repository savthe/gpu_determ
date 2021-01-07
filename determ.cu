/*
 * determ.cu
 *
 *  Created on: 17.02.2010
 *      Author: serge
 */
#include <fstream>
#include <iostream>
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

float* h_f;
float* d_f;
float* h_inverse_integral;
float* d_inverse_integral;
float* h_direct_integral;
float* d_direct_integral;

void init_f(const VelocityGrid& h_vgrid)
{
	h_f = new float[opts.nxyz*opts.npx];
	const float u1 = 1.75;
	const float v1 = 0.25;

	for (int i = 0, index = 0; i < opts.nx; i++) 
		for (int j = 0; j < opts.ny; j++) 
			for (int k = 0; k < opts.nz; k++, index++) 
			{
				float u = h_vgrid.u[i];
				float v = h_vgrid.v[j];
				float w = h_vgrid.w[k];
	
				for (int ix = 0; ix < opts.npx; ix++) 
					h_f[index * opts.npx + ix] = exp (-(u-u1)*(u-u1) - (v-v1)*(v-v1) - w*w); 
			}

	cudaMalloc ((void **) &d_f, opts.nxyz * opts.npx * sizeof (float));
	cudaMemcpy (d_f, h_f, opts.nxyz * opts.npx * sizeof (float), cudaMemcpyHostToDevice);
}

void print_moments(const VelocityGrid& g, const float* f, const float* direct_integral, const float* inverse_integral)
{
  /* EVALUATION OF MOMENTS OF COILLISION INTEGRAL */ 
	float mom0 = 0, momU = 0, momV = 0, momW = 0, mom2 = 0;
	for (int i = 0, index = 0; i < N_X; i++) 
	for (int j = 0; j < N_Y; j++)
	for (int k = 0; k < N_Z; k++, index++) 
	{
		float ci = - f[index] * direct_integral[index] + inverse_integral[index];
		float u = g.u[i];
		float v = g.v[j];
		float w = g.w[k];
		mom0 += ci * g.d3v;
		momU += ci * u * g.d3v;
		momV += ci * v * g.d3v;
		momW += ci * w * g.d3v;
		mom2 += ci * (u*u + v*v + w*w) * g.d3v;
	}

	printf ("M0=%e, MU=%e, MV=%e, MW=%e, M2=%e\n", mom0, momU, momV, momW, mom2);
  /* EVALUATION OF MOMENTS OF COILLISION INTEGRAL */ 
}

void print_results(const VelocityGrid& g, const float* f, const float* direct_integral, const float* inverse_integral)
{
	std::ofstream fresults("results.out");
	for (int index = 0; index < opts.nxyz; index ++) 
		if (g.u_index[index] < g.n_u/2 && g.v_index[index] < g.n_v/2 && g.w_index[index] < g.n_w/2)
			fresults << index << ' ' // index
			<< index / N_YZ << ' ' // I
			<< (index % N_YZ) / N_Z << ' ' // J
	   		<< (index % N_YZ) % N_Z << ' ' // K
			<< f[index] << ' ' // F
			<< direct_integral[index] << ' ' // DC
			<< inverse_integral[index] << ' ' // IC
	   		<< -f[index]*direct_integral[index] + inverse_integral[index] << ' ' // CI
	   		<< -f[index]*direct_integral[index]/inverse_integral[index]+1.0 << std::endl; //ACCURACY
}

void init_integrals()
{
	h_inverse_integral = new float[opts.nxyz*opts.npx];
	h_direct_integral = new float[opts.nxyz*opts.npx];

	cudaMalloc ((void **) &d_direct_integral, opts.nxyz * opts.npx * sizeof (float));
	cudaMalloc ((void **) &d_inverse_integral, opts.nxyz * opts.npx * sizeof (float));
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

	VelocityGrid vgrid(n_points, v_min, v_max, R);

	init_f(vgrid);
	init_integrals();

	float * correction_array;
	init_correction_array (&correction_array, vgrid.device(), opts);

	float * a;
	float * b;
  	init_matrices (vgrid.device(), &b, &a, opts);


	float h_time = 0;
	float * d_time;
	cudaMalloc ((void **) &d_time,  sizeof (float));
	cudaMemcpy (d_time, &h_time, sizeof (float), cudaMemcpyHostToDevice);

  
	GpuTimer timer;
	timer.start();

#if 1  
  
	constexpr int out_step = 5;
	
	for (int step = 0; step <= 0; ++step) 
	{
		cudaMemcpy (&h_time, d_time, sizeof (float), cudaMemcpyDeviceToHost);
		std::cout << "STEP=" << step << ",  TIME=" << h_time << std::endl;
		if (step % out_step == 0) 
			cudaMemcpy (h_f, d_f, opts.nxyz * opts.npx * sizeof (float), cudaMemcpyDeviceToHost);

		evolve_colisions (d_f, d_direct_integral, d_inverse_integral, b, a, correction_array, opts);
		relax_f (d_f, d_direct_integral, d_inverse_integral, 0.5, d_time, opts);
	}
#endif
  
	cudaThreadSynchronize ();

	timer.stop();

	std::cout << "COLLISIONS EVOLVED: " << cudaGetErrorString (cudaGetLastError ()) << std::endl;

	cudaMemcpy (h_direct_integral, d_direct_integral, opts.nxyz * sizeof (float), cudaMemcpyDeviceToHost);
	cudaMemcpy (h_inverse_integral, d_inverse_integral, opts.nxyz * sizeof (float), cudaMemcpyDeviceToHost);


	print_results(vgrid, h_f, h_direct_integral, h_inverse_integral);
	print_moments(vgrid, h_f, h_direct_integral, h_inverse_integral);


  printf ("COLLISION INTEGRAL EVALUATION TOOK %f MS!\n", timer.elapsed());

  cudaFree (d_f);
  cudaFree (d_direct_integral);
  cudaFree (d_inverse_integral);
  cudaFree (correction_array);

  free (h_f);
  free (h_direct_integral);
  free (h_inverse_integral);

  cudaFree (b);
  cudaFree (a);

//  free_device_velocity_grid (d_vgrid);
  //free_host_velocity_grid (h_vgrid);

	printf ("END_PROGRAM: %s\n", 
	cudaGetErrorString (cudaGetLastError ()));
}
