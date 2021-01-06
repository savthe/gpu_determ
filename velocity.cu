/*
 * velocity.cu
 *
 *  Created on: 17.02.2010
 *      Author: serge
 */

#include <assert.h>
#include <cuda.h>
#include <stdio.h>

#include "velocity.h"
#include "options.hpp"

typedef unsigned int uint;

__global__ void fill_correction_array (float * correction_array, const VelocityGrid vgrid, const Options opts)
{
  int i, j, k, index, index1;
  float xi_x, xi_y, xi_z, xi_2;

  for (i = blockIdx.x; i < opts.nx; i += gridDim.x) {
    xi_x = vgrid.u[i];
    for (index1 = threadIdx.x; index1 < N_YZ; index1+= blockDim.x) {
      j = index1 / N_Z;
      k = index1 % N_Z;
      
      xi_y = vgrid.v[j];
      xi_z = vgrid.w[k];
      xi_2 = xi_x * xi_x + xi_y * xi_y + xi_z * xi_z;
      
      index = i * N_YZ + index1;

      XI_X (correction_array)[index] = xi_x;
      XI_Y (correction_array)[index] = xi_y;
      XI_Z (correction_array)[index] = xi_z;
      XI_2 (correction_array)[index] = xi_2;
      XI_X_X (correction_array)[index] = xi_x * xi_x;
      XI_X_Y (correction_array)[index] = xi_x * xi_y;
      XI_X_Z (correction_array)[index] = xi_x * xi_z;
      XI_X_2 (correction_array)[index] = xi_x * xi_2;
      XI_Y_Y (correction_array)[index] = xi_y * xi_y;
      XI_Y_Z (correction_array)[index] = xi_y * xi_z;
      XI_Y_2 (correction_array)[index] = xi_y * xi_2;
      XI_Z_Z (correction_array)[index] = xi_z * xi_z;
      XI_Z_2 (correction_array)[index] = xi_z * xi_2;
      XI_2_2 (correction_array)[index] = xi_2 * xi_2;
    }
  }
}


void init_correction_array (float ** correction_array, const VelocityGrid& vgrid, const Options& opts)
{
  dim3 db, dg;

  cudaMalloc ((void **) correction_array, 
	      opts.nxyz * 14 * sizeof (float));

  db.x = opts.nx;
  db.y = 1;
  db.x = 512;
  db.y = 1;
  db.z = 1;
  fill_correction_array <<<dg, db>>> (*correction_array, vgrid, opts);
}


void init_device_velocity_grid (const VelocityGrid& h_vgrid, VelocityGrid& d_vgrid)
{
	d_vgrid = h_vgrid;


	cudaMalloc ((void **) &d_vgrid.u_index, h_vgrid.n_pnt * sizeof (float));
	cudaMemcpy (d_vgrid.u_index, h_vgrid.u_index, h_vgrid.n_pnt * sizeof(float), cudaMemcpyHostToDevice); 

	cudaMalloc ((void **) &d_vgrid.v_index, h_vgrid.n_pnt * sizeof (float));
	cudaMemcpy (d_vgrid.v_index, h_vgrid.v_index, h_vgrid.n_pnt * sizeof(float), cudaMemcpyHostToDevice); 

	cudaMalloc ((void **) &d_vgrid.w_index, h_vgrid.n_pnt * sizeof (float));
	cudaMemcpy (d_vgrid.w_index, h_vgrid.w_index, h_vgrid.n_pnt * sizeof(float), cudaMemcpyHostToDevice); 

	cudaMalloc ((void **) &d_vgrid.u, h_vgrid.n_u * sizeof (float));
	cudaMemcpy (d_vgrid.u, h_vgrid.u, h_vgrid.n_u * sizeof(float), cudaMemcpyHostToDevice); 

	cudaMalloc ((void **) &d_vgrid.w, h_vgrid.n_w * sizeof (float));
	cudaMemcpy (d_vgrid.w, h_vgrid.w, h_vgrid.n_w * sizeof(float), cudaMemcpyHostToDevice); 

	cudaMalloc ((void **) &d_vgrid.v, h_vgrid.n_v * sizeof (float));
	cudaMemcpy (d_vgrid.v, h_vgrid.v, h_vgrid.n_v * sizeof(float), cudaMemcpyHostToDevice); 
}


void free_host_velocity_grid (VelocityGrid& vgrid)
{
}

void free_device_velocity_grid (VelocityGrid& vgrid)
{
}
