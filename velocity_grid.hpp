#ifndef VELOCITY_GRID_HPP
#define VELOCITY_GRID_HPP

#include "common.hpp"

class VelocityGrid
{
public:
	void init(IntVector n_points, DVector v_min, DVector v_max, float R);
	int n_pnt;
	int* int_params;
//	float* float_params;
	float* w;
	int n_w;
	float* v;
	int n_v;
	float* u;
	int n_u;
//	int* vindex;
	int* u_index;
	int* v_index;
	int* w_index;
	int n_uvw;

	float min_u, min_v, min_w;
	float max_u, max_v, max_w;
	float d_u, d_v, d_w;
	float d3v;
};
/*
typedef struct _VelocityGrid VelocityGrid;
struct _VelocityGrid {
  int n_pnt, n_u, n_v, n_w, n_uvw, n_g;
  int * u_index;
  int * v_index;
  int * w_index;
  int * g_index1;
};

typedef struct _VelocityGridX VelocityGridX;
struct _VelocityGridX {
  short int n_u, n_v, n_w, n_g;
  int aligned_size, aligned_size_g, aligned_vw;
  short int * n_vw;
  short int * vw_index;
  short int * vw_index1;
  short int * start_index_vw;
  short int * ng_vw;
  short int * gindex_vw;
  short int * start_gindex_vw;
  short int * j_range;
  short int * k_range;
};
*/

#endif // VELOCITY_GRID_HPP

