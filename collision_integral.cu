/*
 * collision_integral.cu
 *
 *  Created on: 17.02.2010
 *      Author: serge
 */

//#include <math.h>
//#include <assert.h>
#include <cuda.h>
#include <stdio.h>
#include <assert.h>
#include "velocity.h"
#include "collision_integral.h"
#include <time.h>
#include <device_functions.h>
#include "options.hpp"

__global__ void reset_float_array (float * float_array,
				   int size)
{
  int index;
  
  if (blockIdx.y == 0)
    for (index = threadIdx.x + blockDim.x * blockIdx.x;
	 index < size; index += blockDim.x * gridDim.x)
      float_array[index] = 0.;
}

/*
__global__ void gather_frequency (float * to,
				  float * from)
{
  __shared__ float buf[N_YZ];
  register int i, i1;


  for (i = threadIdx.x; i < N_YZ; i+= blockDim.x)
    buf[i] = 0.0f;
  

  for (i = threadIdx.x; i < N_YZ; i += blockDim.x)
    for (i1 = 0; i1 < N_X; i1++)
      buf[i] += from[(i1 * N_X + blockIdx.x) * N_YZ + i];

  for (i = threadIdx.x; i < N_YZ; i += blockDim.x)
    to[blockIdx.x * N_YZ + i] = buf[i];
  
}
*/


__global__ void evolve_direct_colisions (float * f,
					 float * b,
					 float * direct_integral, 
					 const Options opts)
{
  int index, index_jk;
  short int i, j, k, i1, gi;
  //  __shared__ float f_shared[N_YZ];
  __shared__ float b_shared[N_YZ];
  //  __shared__ float integral_shared[N_YZ];

  i = blockIdx.x;
  index_jk = blockIdx.y;
  j = index_jk / N_Z;
  k = index_jk % N_Z;

  for (i1 = 0; i1 < opts.nx; i1++) {
    short int index1, j1, k1;
    __syncthreads();

#if 0
    for (index = threadIdx.x; index < N_YZ; index+= blockDim.x)
      f_shared[index] = f[i1 * N_YZ + index];
#endif

    gi = (i > i1) ? i - i1 : i1 - i;

    for (index = threadIdx.x; index < N_YZ; index+= blockDim.x)
      b_shared[index] = b[gi * N_YZ + index];

    __syncthreads ();

    for (j1 = 0, index1 = 0; j1 < N_Y; j1++) {
      for (k1 = 0; k1 < N_Z; k1++, index1++) {
	short int gj, gk;

	gj = (j > j1) ? j - j1 : j1 - j;
	gk = (k > k1) ? k - k1 : k1 - k;

	direct_integral[(i * N_YZ + index_jk) * opts.npx + threadIdx.x] += b_shared[gj * N_Z + gk] * 
	  f[(i1 * N_YZ + index1) * opts.npx + threadIdx.x];
      }
    }
#if 0
    direct_integral[(i + i1 * N_X) * N_YZ + index] =
      integral_shared[index];
#endif
  }
}


/*
__global__ void gather_inverse_integral (float * to,
					 float * from)
{
  __shared__ float buf[N_YZ];
  register int i, i1;


  for (i = threadIdx.x; i < N_YZ; i+= blockDim.x)
    buf[i] = 0.0f;
  

  for (i = threadIdx.x; i < N_YZ; i += blockDim.x)
    for (i1 = 0; i1 < N_X * (N_X >> 1); i1++)
      buf[i] += from[(i1 * N_X + blockIdx.x) * N_YZ + i];

  for (i = threadIdx.x; i < N_YZ; i += blockDim.x)
    to[blockIdx.x * N_YZ + i] = buf[i];
}
*/

__global__ void evolve_inverse_colisions (float * f,
					  float * a,
					  float * inverse_integral,
					  short int in1,
					  short int jn1,
					  short int kn1,
					  const Options opts)
{
  short int i, i1, gi, in, il, il1;
  __shared__ float a_shared[N_YZ];
  //  __shared__ float integral_shared[opts.npx];
  short int index, index1;

  //  i1 = blockIdx.y;
  i = blockIdx.x;
  index1 = blockIdx.y;

  for (i1 = 0; i1 < opts.nx; i1++) {
    in = (in1 << 1) | ((i + i1) & 1);

    il = (i + i1 + in) >> 1;
    if (il > (opts.nx - 1)) continue;
    il1 = i + i1 - in;
    if (il1 < 0) continue;
    il1 = (i + i1 - in ) >> 1;


    gi = (i > i1) ? (i - i1) : (i1 - i);

    index = jn1 * (N_Z >> 1) + kn1;

    short int l_index, l1_index;                            

    for (l_index = threadIdx.x; l_index < N_YZ; l_index += blockDim.x) // copy array a to shared 
      a_shared[l_index] = a[((in1 * opts.nx + gi) * (N_YZ >> 2) + index) 
			   * N_YZ + l_index]; // memory
    __syncthreads ();
  
    register short int j, k, j1, k1, gj, gk, jl, jl1, kl, kl1, jn, kn;

    j = index1 / N_Z; // defines the point i in which collision integral is computed
    k = index1 % N_Z;
#if 1	
    for (j1 = 0; j1 < N_Y; j1 ++) {//the sum n_ + i_ + j_ must be even
      gj = (j > j1) ? (j - j1) : (j1 - j);
      jn = (gj & 1) | (jn1 << 1);
      jl = (j + j1 + jn) >> 1; // the indices 
      jl1 = j + j1 - jn;  // of post-collision

      for (k1 = 0; k1 < N_Z; k1 ++) { // summation on j loop
	gk = (k > k1) ? (k - k1) : (k1 - k);
	kn = (gk & 1) | (kn1 << 1);
	
	int a_index = gj * N_Z + gk;       //get 1 index from 2 components
	
#if 1	
	kl = (k + k1 + kn) >> 1; // velocities
	kl1 = k + k1 - kn;
#endif
        
	l_index = (jl < N_Y  && kl < N_Z) ? jl * N_Z + kl : N_YZ; 
	l1_index = (jl1 >= 0  && kl1 >= 0) ? (jl1 >> 1) * N_Z + (kl1 >> 1) : N_YZ; 
#if 0
	integral_shared[threadIdx.x] = 
	  (l_index == N_YZ || l1_index == N_YZ) ? 0.0f :
	  a_shared[a_index] * 
	  (f[(il * N_YZ + l_index)*opts.npx + threadIdx.x] * f[(il1 * N_YZ + l1_index)*opts.npx + threadIdx.x] + // accumulating the collision 
	   f[(il * N_YZ + l1_index)*opts.npx + threadIdx.x] * f[(il1 * N_YZ + l_index)*opts.npx + threadIdx.x]); // integral
	
	l_index = (jl < N_Y  && kl1 >= 0) ? jl * N_Z + (kl1 >> 1) : N_YZ; 
	l1_index = (jl1 >= 0  && kl < N_Z) ? (jl1 >> 1) * N_Z + kl : N_YZ;

	integral_shared[threadIdx.x] += 
	  (l_index == N_YZ || l1_index == N_YZ) ? 0.0f :
	  a_shared[a_index] * 
	  (f[(il * N_YZ + l_index)*opts.npx + threadIdx.x] * f[(il1 * N_YZ + l1_index)*opts.npx + threadIdx.x] + // accumulating the collision 
	   f[(il * N_YZ + l1_index)*opts.npx + threadIdx.x] * f[(il1 * N_YZ + l_index)*opts.npx + threadIdx.x]); // integral

	inverse_integral[(i * N_X + index1) * opts.npx + threadIdx.x] += integral_shared[threadIdx.x];
#elif 1
	float tmp = 
	  (l_index == N_YZ || l1_index == N_YZ) ? 0.0f :
	  a_shared[a_index] * 
	  (f[(il * N_YZ + l_index)*opts.npx + threadIdx.x] * f[(il1 * N_YZ + l1_index)*opts.npx + threadIdx.x] + // accumulating the collision 
	   f[(il * N_YZ + l1_index)*opts.npx + threadIdx.x] * f[(il1 * N_YZ + l_index)*opts.npx + threadIdx.x]); // integral
	
	l_index = (jl < N_Y  && kl1 >= 0) ? jl * N_Z + (kl1 >> 1) : N_YZ; 
	l1_index = (jl1 >= 0  && kl < N_Z) ? (jl1 >> 1) * N_Z + kl : N_YZ;

	tmp += 
	  (l_index == N_YZ || l1_index == N_YZ) ? 0.0f :
	  a_shared[a_index] * 
	  (f[(il * N_YZ + l_index)*opts.npx + threadIdx.x] * f[(il1 * N_YZ + l1_index)*opts.npx + threadIdx.x] + // accumulating the collision 
	   f[(il * N_YZ + l1_index)*opts.npx + threadIdx.x] * f[(il1 * N_YZ + l_index)*opts.npx + threadIdx.x]); // integral

	inverse_integral[(i * opts.nx + index1) * opts.npx + threadIdx.x] += tmp;
#elif 0
	inverse_integral[(i * opts.nx + index1) * opts.npx + threadIdx.x] += 
	  (l_index == N_YZ || l1_index == N_YZ) ? 0.0f :
	  a_shared[a_index] * 
	  (f[(il * N_YZ + l_index)*opts.npx + threadIdx.x] * f[(il1 * N_YZ + l1_index)*opts.npx + threadIdx.x] + // accumulating the collision 
	   f[(il * N_YZ + l1_index)*opts.npx + threadIdx.x] * f[(il1 * N_YZ + l_index)*opts.npx + threadIdx.x]); // integral
	
	l_index = (jl < N_Y  && kl1 >= 0) ? jl * N_Z + (kl1 >> 1) : N_YZ; 
	l1_index = (jl1 >= 0  && kl < N_Z) ? (jl1 >> 1) * N_Z + kl : N_YZ;

	inverse_integral[(i * opts.nx + index1) * opts.npx + threadIdx.x] += 
	  (l_index == N_YZ || l1_index == N_YZ) ? 0.0f :
	  a_shared[a_index] * 
	  (f[(il * N_YZ + l_index)*opts.npx + threadIdx.x] * f[(il1 * N_YZ + l1_index)*opts.npx + threadIdx.x] + // accumulating the collision 
	   f[(il * N_YZ + l1_index)*opts.npx + threadIdx.x] * f[(il1 * N_YZ + l_index)*opts.npx + threadIdx.x]); // integral
#endif
      }
    }
#endif
  }
#if 0
  for (index = threadIdx.x; index < N_YZ; index+= blockDim.x) 
    inverse_integral[(i + i1 * N_X) * N_YZ + index] += 
      integral_shared[index];  // copy the results from shared memory to global memory  
#endif
}

/*
void check_a_array (float * a)
{
  float * h_a;
  short int in1, jn1, kn1, gi, gj, gk;
  int index;
  h_a = (float *) malloc (opts.nxyz * (opts.nxyz >> 3) * sizeof (float));
  cudaMemcpy (h_a, a, 
	      opts.nxyz * (opts.nxyz >> 3) * sizeof (float),
	      cudaMemcpyDeviceToHost);
  for (in1 = 0, index = 0; in1 < (N_X >> 1); in1++)
    for (gi = 0; gi < N_X; gi++)
      for (jn1 = 0; jn1 < (N_Y >>1); jn1++)
	for (kn1 = 0; kn1 < (N_Z >>1); kn1++)
	  for (gj = 0; gj < N_Y; gj++)
	    for (gk = 0; gk < N_X; gk++, index++)
	      //	      if (! (h_a[index] >= 0 && (h_a[index] < 2)))
	      if (gi < 3 && gj < 3 && gk < 3 && in1 == 0 && jn1 == 0 && 
		  kn1 == 1 && h_a[index] > ERR)
		printf ("INDICES:%d %d %d %d %d %d %d A=%e\n",
			index, in1, gi, jn1, kn1, gj, gk, 
			h_a[index]);
  free(h_a);
} 
*/


__global__ void before_correction(float * inverse_integral,
				  float * frequency,
				  float * f,
				  float * correction_array,
	//			  float A[25*NP_X],
				  float * A,
				  float B[5*NP_X],
				  const Options opts)
{
  __shared__ float shared_buf[1024];
  int i, index;
  int ix = blockIdx.y;

  for (index = threadIdx.x; index < 1024; index += blockDim.x)
    shared_buf[index] = 0;

  __syncthreads ();

  if (blockIdx.x == 0) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += inverse_integral[index + ix*opts.nxyz] - 
	f[index + ix*opts.nxyz] * frequency[index + ix*opts.nxyz];
  }
#if 1
  else if (blockIdx.x == 1) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += (inverse_integral[index + ix*opts.nxyz] - 
			    f[index + ix*opts.nxyz] * frequency[index + ix*opts.nxyz]) * 
	(XI_X (correction_array)[index]);
  }
  else if (blockIdx.x == 2) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += (inverse_integral[index + ix*opts.nxyz] -
			    f[index + ix*opts.nxyz] * frequency[index + ix*opts.nxyz]) * 
	(XI_Y (correction_array)[index]);
  }
  else if (blockIdx.x == 3) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += (inverse_integral[index + ix*opts.nxyz] -
			    f[index + ix*opts.nxyz] * frequency[index + ix*opts.nxyz]) * 
	(XI_Z (correction_array)[index]);
  }
  else if (blockIdx.x == 4) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += (inverse_integral[index + ix*opts.nxyz] - 
			    f[index + ix*opts.nxyz] * frequency[index + ix*opts.nxyz]) * 
	(XI_2 (correction_array)[index]);
  }
  else if (blockIdx.x == 5) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz];
  }
  else if (blockIdx.x == 6) {
#if 1
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz] *
	(XI_X (correction_array)[index]);
#else
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += 
	(XI_X (correction_array)[index]);
#endif
  }
  else if (blockIdx.x == 7) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz] *
	(XI_Y (correction_array)[index]);
  }
  else if (blockIdx.x == 8) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz] *
	(XI_Z (correction_array)[index]);
  }
  else if (blockIdx.x == 9) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz] *
	(XI_2 (correction_array)[index]);
  }
  else if (blockIdx.x == 10) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz] *
	(XI_X_X (correction_array)[index]);
  }
  else if (blockIdx.x == 11) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz] *
	(XI_X_Y (correction_array)[index]);
  }
  else if (blockIdx.x == 12) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz] *
	(XI_X_Z (correction_array)[index]);
  }
  else if (blockIdx.x == 13) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz] *
	(XI_X_2 (correction_array)[index]);
  }
  else if (blockIdx.x == 14) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz] *
	(XI_Y_Y (correction_array)[index]);
  }
  else if (blockIdx.x == 15) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz] *
	(XI_Y_Z (correction_array)[index]);
  }
  else if (blockIdx.x == 16) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz] *
	(XI_Y_2 (correction_array)[index]);
  }
  else if (blockIdx.x == 17) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz] *
	(XI_Z_Z (correction_array)[index]);
  }
  else if (blockIdx.x == 18) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz] *
	(XI_Z_2 (correction_array)[index]);
  }
  else if (blockIdx.x == 19) {
    for (index = threadIdx.x; index < opts.nxyz; index += blockDim.x)
#if 1
      shared_buf[index & 1023] += frequency[index + ix*opts.nxyz] * f[index + ix*opts.nxyz] *
	(XI_2_2 (correction_array)[index]);
#else
      shared_buf[index & 1023] += 
	(XI_2_2 (correction_array)[index]);
#endif 
  }
  __syncthreads ();

#if 1
  i = 512;
  while (i > 0) {
    for (index = threadIdx.x; index < i; index += blockDim.x)
      shared_buf [index] += shared_buf [index + i];
    i = i >> 1;
    __syncthreads ();
  }
#endif

  if (blockIdx.x == 0){
    if (threadIdx.x == 0) B[0 + 5*ix] = shared_buf[0];
  }
  else if (blockIdx.x == 1) {
    if (threadIdx.x == 0) B[1 + 5*ix] = shared_buf[0];
  }
  else if (blockIdx.x == 2) {
    if (threadIdx.x == 0) B[2 + 5*ix] = shared_buf[0];
  }
  else if (blockIdx.x == 3) {
    if (threadIdx.x == 0) B[3 + 5*ix] = shared_buf[0];
  }
  else if (blockIdx.x == 4) {
    if (threadIdx.x == 0) B[4 + 5*ix] = shared_buf[0];
  }
  else if (blockIdx.x == 5) {
    if (threadIdx.x == 0) A[0 + 25*ix] = shared_buf[0];
  }
  else if (blockIdx.x == 6) {
    if (threadIdx.x == 0) {
      A[1 * 5 + 0 + 25*ix] = shared_buf[0];
      A[0 * 5 + 1 + 25*ix] = shared_buf[0];
    }
  }
  else if (blockIdx.x == 7) {
    if (threadIdx.x == 0) {
      A[2 * 5 + 0 + 25*ix] = shared_buf[0];
      A[0 * 5 + 2 + 25*ix] = shared_buf[0];
    }
  }
  else if (blockIdx.x == 8) {
    if (threadIdx.x == 0) {
      A[3 * 5 + 0 + 25*ix] = shared_buf[0];
      A[0 * 5 + 3 + 25*ix] = shared_buf[0];
    }
  }
  else if (blockIdx.x == 9) {
    if (threadIdx.x == 0) {
      A[4 * 5 + 0 + 25*ix] = shared_buf[0];
      A[0 * 5 + 4 + 25*ix] = shared_buf[0];
    }
  }
  else if (blockIdx.x == 10) {
    if (threadIdx.x == 0) {
      A[1 * 5 + 1 + 25*ix] = shared_buf[0];
    }
  }
  else if (blockIdx.x == 11) {
    if (threadIdx.x == 0) {
      A[1 * 5 + 2 + 25*ix] = shared_buf[0];
      A[2 * 5 + 1 + 25*ix] = shared_buf[0];
    }
  }
  else if (blockIdx.x == 12) {
    if (threadIdx.x == 0) {
      A[1 * 5 + 3 + 25*ix] = shared_buf[0];
      A[3 * 5 + 1 + 25*ix] = shared_buf[0];
    }
  }
  else if (blockIdx.x == 13) {
    if (threadIdx.x == 0) {
      A[1 * 5 + 4 + 25*ix] = shared_buf[0];
      A[4 * 5 + 1 + 25*ix] = shared_buf[0];
    }
  }
  else if (blockIdx.x == 14) {
    if (threadIdx.x == 0) {
      A[2 * 5 + 2 + 25*ix] = shared_buf[0];
    }
  }
  else if (blockIdx.x == 15) {
    if (threadIdx.x == 0) {
      A[2 * 5 + 3 + 25*ix] = shared_buf[0];
      A[3 * 5 + 2 + 25*ix] = shared_buf[0];
    }
  }
  else if (blockIdx.x == 16) {
    if (threadIdx.x == 0) {
      A[2 * 5 + 4 + 25*ix] = shared_buf[0];
      A[4 * 5 + 2 + 25*ix] = shared_buf[0];
    }
  }
  else if (blockIdx.x == 17) {
    if (threadIdx.x == 0) {
      A[3 * 5 + 3 + 25*ix] = shared_buf[0];
    }
  }
  else if (blockIdx.x == 18) {
    if (threadIdx.x == 0) {
      A[3 * 5 + 4 + 25*ix] = shared_buf[0];
      A[4 * 5 + 3 + 25*ix] = shared_buf[0];
    }
  }
  else if (blockIdx.x == 19) {
    if (threadIdx.x == 0) {
      A[4 * 5 + 4 + 25*ix] = shared_buf[0];
    }
  }
#endif
}


__global__ void linear_solver(const float A[25*NP_X], 
			      const float B[5*NP_X], float X[5*NP_X])  
{
  __shared__ float A1[25], B1[5];
  int i, j, k, ix;
  ix = blockIdx.x;

  if (blockIdx.x > 0) return;
  
  i = threadIdx.x;
  j = threadIdx.y;
  if (j < 5 && i < 5) {
    A1[i * 5 + j] = A[i * 5 + j + ix*5];
    if (j == 0)
	B1[i] = B[i + ix*5];
  }
  __syncthreads ();
  for (k = 0; k < 4; k++) {
    if (i < 5 && j < 5){
      if (i > k){
	if (j >= k)
	  A1[i * 5 + j] -= A1[k * 5 + j] * A1[i * 5 +k] / A1[k * 5 + k];
	}
      else if (i == k) {
	if (j > k)
	  B1[j] -= B1[k] * A1[j * 5 + k] / A1[k * 5 + k];
      }
    }
      __syncthreads ();
  }

  for (k = 4; k > 0; k--) {
    if (i < k && j == 0)
      B1[i] -= B1[k] * A1[i * 5 + k] / A1[k * 5 + k];    
    __syncthreads ();
  }

  if (i < 5 && j == 0)
    X[i + ix*5] = B1[i] / A1[i * 5 + i];
}


__global__ void correct_frequency (float * direct_integral,
				   const float * correction_array,
				   float X[5*NP_X],
				   const Options opts)
{
  __shared__ float pl_shared[512];
  int index;
  int ix = blockIdx.y;

  if (blockIdx.x >= opts.nx) return;

  for (index = threadIdx.x; index < N_YZ; index += blockDim.x) {
    int index1 = blockIdx.x * N_YZ + index;
    pl_shared[threadIdx.x] = 1.0 + X[0 + 5*ix];
    pl_shared[threadIdx.x] += X[1 + 5*ix] * XI_X (correction_array) [index1];
    pl_shared[threadIdx.x] += X[2 + 5*ix] * XI_Y (correction_array) [index1];
    pl_shared[threadIdx.x] += X[3 + 5*ix] * XI_Z (correction_array) [index1];
    pl_shared[threadIdx.x] += X[4 + 5*ix] * XI_2 (correction_array) [index1];
    direct_integral[index1 + ix*opts.nxyz] *= pl_shared[threadIdx.x];
  }
}

__global__ void transpose_float (float * from, float * to, int rows, int cols)
{
  __shared__ float A[16*17];

  int i, j;

  j = threadIdx.x + blockIdx.x * blockDim.x;
  i = threadIdx.y + blockIdx.y * blockDim.y;

  if (i < rows && j < cols)
    A[threadIdx.x + threadIdx.y * 17] = from [j + i * cols];

  __syncthreads();

  i = threadIdx.x + blockIdx.y * blockDim.y;
  j = threadIdx.y + blockIdx.x * blockDim.x;
  
  if (i < rows && j < cols)
    to[i + j * rows] = A[threadIdx.y + 17 * threadIdx.x]; 
}


void evolve_colisions (float * f,
		       float * direct_integral,
		       float * inverse_integral,
		       float * b,
		       float * a,
		       float * correction_array,
			   const Options& opts)
{
  cudaEvent_t start, stop;
  float elapsedTime;
  dim3 dg1;
  dim3 db1;
  short int in1, jn1, kn1;
  //  float * integral1;
  float * X;
  float * A;
  float * B;

  static float * ft = NULL;
  static float *direct_t = NULL;
  static float *inverse_t = NULL;
  
  if (ft == NULL)
    cudaMalloc ((void **) &ft, opts.npx * opts.nxyz * sizeof (float));
  if (direct_t == NULL)
    cudaMalloc ((void **) &direct_t, opts.npx * opts.nxyz * sizeof (float));
  if (inverse_t == NULL)
    cudaMalloc ((void **) &inverse_t, opts.npx * opts.nxyz * sizeof (float));


  float A1[opts.npx * 25], B1[opts.npx * 5], X1[opts.npx * 5];
  
#if 0
  check_a_array(a);
#endif


  dg1.x = opts.nx;
  dg1.y = N_YZ;
  db1.x = opts.npx;
  db1.y = 1;
  db1.z = 1;

  reset_float_array <<<1024, 512>>> (direct_integral, opts.npx * opts.nxyz);
  reset_float_array <<<1024, 512>>> (direct_t, opts.npx * opts.nxyz);
  reset_float_array <<<1024, 512>>> (inverse_t, opts.npx * opts.nxyz);
  reset_float_array <<<1024, 512>>> (ft, opts.npx * opts.nxyz);
  cudaThreadSynchronize ();
  printf ("RESET DONE!, %s\n", 
	  cudaGetErrorString (cudaGetLastError ()));

  cudaEventCreate(&start);
  cudaEventCreate(&stop);
#if 1

  cudaEventRecord(start, 0);
  evolve_direct_colisions <<<dg1, db1 >>> (f, b, direct_integral,opts);
  
  cudaThreadSynchronize ();

#if 0
  dg1.y = 1;
  gather_frequency <<<dg1, db1>>>(direct_integral, integral1);
#endif

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  cudaEventElapsedTime(&elapsedTime, start, stop);
  printf ("NEW OPTIMIZED DIRECT COLISIONS COMPUTING TOOK %f MS!, %s\n", 
	  elapsedTime, cudaGetErrorString (cudaGetLastError ()));

  cudaThreadSynchronize();

  printf ("BEFORE INVERSE COLISIONS COMPUTING!, %s\n", 
	  cudaGetErrorString (cudaGetLastError ()));

  reset_float_array <<<1024, 512>>> (inverse_integral, opts.npx * opts.nxyz);
  cudaThreadSynchronize ();
  cudaEventRecord(start, 0);
  dg1.x = opts.nx;
  dg1.y = N_YZ;
  db1.x = opts.npx;
  db1.y = 1;
  time_t now = time(0);
  for (in1 = 0; in1 < (opts.nx >> 1); in1++)
    for (jn1 = 0; jn1 < (N_Y >> 1); jn1++)
      for (kn1 = 0; kn1 < (N_Z >> 1); kn1++) {
	cudaEventRecord(start, 0);
	evolve_inverse_colisions <<<dg1, db1>>> (f, a, inverse_integral, in1, jn1, kn1, opts);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);

	cudaEventElapsedTime(&elapsedTime, start, stop);
	printf ("CALL WITH I=%d, J=%d, K=%d TOOK %f MS!, %s\n", in1, jn1, kn1, 
		elapsedTime, cudaGetErrorString (cudaGetLastError ()));
      }
  printf("SUMMARY KERNEL TIME: %d", time(0) - now);
#endif

  cudaThreadSynchronize ();

#if 0
  dg1.x = N_X;
  dg1.y = 1;
  //  db1.x = 256;
  //  db1.y = 1;

  gather_frequency <<<dg1, db1>>>(inverse_integral, integral1);
#endif

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  printf ("OPTIMIZED INVERSE COLISIONS COMPUTING TOOK %f MS!, %s\n", 
	  elapsedTime, cudaGetErrorString (cudaGetLastError ()));

  cudaMalloc ((void **) &A, 25 * opts.npx * sizeof (float));
  cudaMalloc ((void **) &B, 5 * opts.npx * sizeof (float));
  cudaMalloc ((void **) &X, 5 * opts.npx * sizeof (float));

  dg1.y = opts.nxyz / 16 + 1;
  dg1.x = opts.npx / 16 + 1;
  db1.x = 16;
  db1.y = 16;

#if 1
  transpose_float <<<dg1, db1>>> (f, ft, opts.nxyz, opts.npx);
  cudaThreadSynchronize ();
  printf ("TRANSPOSE1 DONE!, %s\n", 
	  cudaGetErrorString (cudaGetLastError ()));
  transpose_float <<<dg1, db1>>> (inverse_integral, inverse_t, opts.nxyz, opts.npx);
  cudaThreadSynchronize ();
  printf ("TRANSPOSE2 DONE!, %s\n", 
	  cudaGetErrorString (cudaGetLastError ()));
  transpose_float <<<dg1, db1>>> (direct_integral, direct_t, opts.nxyz, opts.npx);
  cudaThreadSynchronize ();
  printf ("TRANSPOSE3 DONE!, %s\n", 
	  cudaGetErrorString (cudaGetLastError ()));
#endif

  dg1.x = 20;
  dg1.y = opts.npx;
  db1.x = 512;
  db1.y = 1;
  db1.z = 1;
  before_correction <<<dg1, db1>>>(inverse_t, direct_t,
		     ft, correction_array, A, B, opts);

  cudaThreadSynchronize ();
  printf ("BEFORE_CORRECTION DONE!, %s\n", 
	  cudaGetErrorString (cudaGetLastError ()));
  cudaMemcpy (A1, A, 25 * opts.npx * sizeof (float),
	      cudaMemcpyDeviceToHost);
  cudaMemcpy (B1, B, 5 * opts.npx * sizeof (float),
	      cudaMemcpyDeviceToHost);

  for (in1 = 0, kn1 = 0; in1 < 5 ; in1++){
    for (jn1 = 0; jn1 < 5; jn1++, kn1++)
      printf ("%e ", A1[kn1]);
    printf("|| %e\n", B1[in1]);
  }

  dg1.x = opts.npx;
  dg1.y = 1;
  db1.x = 5;
  db1.y = 5;
  db1.z = 1;
  linear_solver <<<dg1, db1>>>(A, B, X);

#if 1
  cudaThreadSynchronize ();
  printf ("LINEAR DONE!, %s\n", 
	  cudaGetErrorString (cudaGetLastError ()));

  cudaMemcpy (X1, X, 5 * opts.npx * sizeof (float),
	      cudaMemcpyDeviceToHost);

  for (in1 = 0; in1 < 15 ; in1++){
    float res = 0.0;
    for (jn1 = 0; jn1 < 5; jn1++)
      res += A1[5 * in1 + jn1] * X1[jn1];
    printf ("%e : %e", X1[in1], res);
    printf("|| %e\n", B1[in1]);
  }
#endif

  dg1.x = opts.nx;
  dg1.y = opts.npx;
  db1.x = 512;
  db1.y = 1;
  db1.z = 1;
  correct_frequency <<<dg1, db1>>>(direct_t, 
				   correction_array,
				   X, 
				   opts);
  cudaThreadSynchronize ();
  printf ("CORRECTION DONE!, %s\n", 
	  cudaGetErrorString (cudaGetLastError ()));

  float * hf;
  int i, j, index;

  hf = (float *) malloc (opts.nxyz * opts.npx * sizeof (float));
  cudaMemcpy (hf, direct_t, opts.nxyz * opts.npx * sizeof (float),
	      cudaMemcpyDeviceToHost);

  for (i = 0, index = 0; i < N_YZ; i++) {
    for (j = 0; j < opts.nx; j++, index++) 
      printf ("%7.1e ", hf[index]);
    printf ("\n");  
  }  

  cudaFree (A);
  cudaFree (B);
  cudaFree (X);

  dg1.x = opts.nxyz / 16 + 1;
  dg1.y = opts.npx / 16 + 1;
  db1.x = 16;
  db1.y = 16;

  transpose_float <<<dg1, db1>>> (direct_t, direct_integral, opts.npx, opts.nxyz);
  //  assert(!"here");
}


__global__ void get_delta_t (float * frequency,
			     float courant,
			     float * delta_t,
				 const Options opts)
{
  __shared__ float delta_t_shared[1024];
  int index, i;
  if (blockIdx.x > 0) return;

  for (index = threadIdx.x; index < 1024; index += blockDim.x)
    delta_t_shared[index] = 1.0e10;

  for (index = threadIdx.x; index < opts.nxyz * opts.npx; index += blockDim.x)
    delta_t_shared[index & 1023] =   courant/frequency[index] <
      delta_t_shared[index & 1023] ? courant/frequency[index] :
      delta_t_shared[index & 1023];

  __syncthreads();


  i = 512;
  while (i > 0) {
    for (index = threadIdx.x; index < i; index += blockDim.x)
      delta_t_shared [index] = 
	delta_t_shared[index] < delta_t_shared[index + i] ? 
	delta_t_shared[index] : delta_t_shared[index + i];
    i = i >> 1;
    __syncthreads ();
  }
  if (threadIdx.x == 0)
    *delta_t = delta_t_shared[0];
}

__global__ void evolve_f (float * f,
			  float * frequency,
			  float * inverse_integral,
			  float * delta_t,
			  float * time)
{
  int index, index1;
  __shared__ float df_shared[512];
  __shared__ float dt;

  if (threadIdx.x == 0) dt = (*delta_t);

  __syncthreads ();

  for (index1 = threadIdx.x; index1 < N_YZ; index1 += blockDim.x) {
    index = blockIdx.x * N_YZ + index1;
    df_shared[threadIdx.x] =  f[index] *
      (1.0f - frequency[index] * dt) +
      inverse_integral[index] * dt;
    f[index] = df_shared[threadIdx.x];
  }
  if (blockIdx.x == 0 && threadIdx.x == 0)
    (*time) += dt;
  __syncthreads ();
}

void relax_f (float * f,
	      float * frequency,
	      float * inverse_integral,
	      float courant,
	      float * time,
		  const Options& opts)
{
  float * delta_t;
  dim3 dg, db;
  float dt;
  dg.x = 1;
  dg.y = 1;
  db.x = 512;
  db.y = 1;
  db.z = 1;



  cudaMalloc ((void **) &delta_t, 
	      sizeof (float));


  get_delta_t <<<dg, db>>> (frequency, courant, delta_t, opts);
  cudaThreadSynchronize ();
  cudaMemcpy (&dt, delta_t, 
	      sizeof (float),
	      cudaMemcpyDeviceToHost);
  printf ("AFTER DELTA_T COMPUTING, DELTA_T=%e!, %s\n", dt,
	  cudaGetErrorString (cudaGetLastError ()));


  dg.x = opts.nx*opts.npx;
#if 1
  evolve_f <<<dg, db>>> (f, frequency, inverse_integral, delta_t, 
			 time);
#endif

  cudaMemcpy (&dt, delta_t, 
	      sizeof (float),
	      cudaMemcpyDeviceToHost);

  cudaThreadSynchronize ();
  cudaFree (delta_t);
  printf ("AFTER f COMPUTING, DELTA_T=%e!, %s\n", dt, 
	  cudaGetErrorString (cudaGetLastError ()));
}


