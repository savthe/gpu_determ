#include <cassert>
#include "velocity_grid.hpp"
#include "velocity.h"


bool check_if_point_is_in_domain (IntVector n_points, 
			                     DVector v_min, 
			                     DVector v_max,
			                     int i, int j, int k,
			                     float R)
{
  float u = (v_max.x - v_min.x) * (2 * i + 1 - n_points.x) / 2. / 
              n_points.x; 
  float v = (v_max.y - v_min.y) * (2 * j + 1 - n_points.y) / 2. / 
              n_points.y; 
  float w = (v_max.z - v_min.z) * (2 * k + 1 - n_points.z) / 2. / 
              n_points.z; 
  return (R <= 0 || u * u + v * v + w * w < R * R);
}


void VelocityGrid::init(IntVector n_points, DVector v_min, DVector v_max, float R)
{
	int** int_arr = &int_params;
//	float** d_arr = &float_params;
  int i, j, k, index, index1;
  n_pnt = 0;
  
  assert (n_points.x > 0 && n_points.y > 0 && n_points.z > 0);
  assert (v_min.x < v_max.x && v_min.y < v_max.y && v_min.z < v_max.z);
 
  for (i = 0; i < n_points.x; i++)
    for (j = 0; j < n_points.y; j++)
      for (k = 0; k < n_points.z; k++)
	if (check_if_point_is_in_domain (n_points, v_min, v_max,
					 i, j, k, R))
	  n_pnt++;
  
  (*int_arr) = (int *) malloc ((5 + n_pnt * 4 + 
                                n_points.x * n_points.y * n_points.z) 
                               * sizeof (int));

  N_PNT (*int_arr) = n_pnt;
  //N_U (*int_arr)   = n_points.x;
  n_u = n_points.x;
  //N_V (*int_arr)   = n_points.y;
  n_v = n_points.y;
  //N_W (*int_arr)   = n_points.z;
  n_w = n_points.z;
  N_UVW (*int_arr) = n_points.x * n_points.y * n_points.z;

  
  /*
  (*d_arr) = (float *) malloc ((10 + n_points.x + n_points.y +
				 n_points.z + 
				 n_pnt * 3) * sizeof (float));
				 */
  u = (float*) malloc(n_u*sizeof(float));
  v = (float*) malloc(n_v*sizeof(float));
  w = (float*) malloc(n_points.z*sizeof(float));

  u_index = (int*) malloc(n_pnt*sizeof(int));
  v_index = (int*) malloc(n_pnt*sizeof(int));
  w_index = (int*) malloc(n_pnt*sizeof(int));

//  vindex = (int*) malloc(n_points.x*n_points.y*n_points.z*sizeof(int));

  /*
  MIN_U (*d_arr) = v_min.x;
  MIN_V (*d_arr) = v_min.y;
  MIN_W (*d_arr) = v_min.z;
  MAX_U (*d_arr) = v_max.x;
  MAX_V (*d_arr) = v_max.y;
  MAX_W (*d_arr) = v_max.z;
  D_U (*d_arr) = (MAX_U (*d_arr) - MIN_U (*d_arr)) / n_u;//N_U (*int_arr);
  D_V (*d_arr) = (MAX_V (*d_arr) - MIN_V (*d_arr)) / n_v;//N_V (*int_arr);
  D_W (*d_arr) = (MAX_W (*d_arr) - MIN_W (*d_arr)) / n_w;//N_W (*int_arr);
  D3V (*d_arr) = D_U (*d_arr) * D_V (*d_arr) * D_W(*d_arr);
  */
  min_u = v_min.x;
  min_v = v_min.y;
  min_w = v_min.z;
  max_u = v_max.x;
  max_v = v_max.y;
  max_w = v_max.z;
  d_u = (max_u - min_u) / n_u;//N_U (*int_arr);
  d_v = (max_v - min_v) / n_v;//N_V (*int_arr);
  d_w = (max_w - min_w) / n_w;//N_W (*int_arr);
  d3v = d_u*d_v*d_w;
  
  
  for (i = 0; i < n_u; i++)
    u[i] = min_u + d_u * (2.*i+1.)/2.;
  for (i = 0; i < n_v; i++)
    v[i] = min_v + d_v * (2.*i+1.)/2.;
    //V_1 (*d_arr, *int_arr)[i] = MIN_V (*d_arr) + D_V (*d_arr) * (2.*i+1.)/2.;
  for (i = 0; i < n_w; i++)
    w[i] = min_w + d_w * (2.*i+1.)/2.;
    //W_1 (*d_arr, *int_arr)[i] = MIN_W (*d_arr) + D_W (*d_arr) * (2.*i+1.)/2.;
  
  index = 0;
  for (i = 0, index1 = 0; i < n_u; i++)
    for (j = 0; j < n_v; j++)
      for (k = 0; k < n_w; k++, index1++)
	    if (check_if_point_is_in_domain (n_points, v_min, v_max,
	                                      i, j, k, R)) {
	      //U_0 (*d_arr, *int_arr)[index] = U_1 (*d_arr, *int_arr)[i];
	      //V_0 (*d_arr, *int_arr)[index] = V_1 (*d_arr, *int_arr)[j];
	      //W_0 (*d_arr, *int_arr)[index] = W_1 (*d_arr, *int_arr)[k];
		  //w0[index] = W_1(*d_arr, *int_arr)[k];
	      u_index[index] = i;
	      v_index[index] = j;
	      //W_INDEX (*int_arr)[index] = k;
		  w_index[index] = k;
		  //vindex[index1] = index;
	      //INDEX (*int_arr)[index1] = index;
//	      INDEX1 (*int_arr)[index] = index1;
	      index++;
	    }
//	    else 
	      //vindex[index1] = -1;
	      //INDEX (*int_arr)[index1] = -1;
  
}

