#ifndef VELOCITY_GRID_HPP
#define VELOCITY_GRID_HPP

#include "common.hpp"

class VelocityGridBase
{
	//VelocityGridBase(const VelocityGridBase&) = delete;
	//VelocityGridBase& operator=(const VelocityGridBase&) = delete;

public:
	VelocityGridBase(){}
	virtual ~VelocityGridBase() {};

	int n_pnt = 0;

	int n_u = 0;
	int n_w = 0;
	int n_v = 0;

	float* u = 0;
	float* v = 0;
	float* w = 0;

	int* u_index = 0;
	int* v_index = 0;
	int* w_index = 0;

	int n_uvw = 0;

	float min_u, min_v, min_w;
	float max_u, max_v, max_w;
	float d_u, d_v, d_w;
	float d3v;
};

class DeviceVelocityGrid;

class VelocityGrid: public VelocityGridBase
{
public:
	VelocityGrid() {}
	VelocityGrid(Vector3i n_points, Vector3f v_min, Vector3f v_max, float r);
	void init(Vector3i n_points, Vector3f v_min, Vector3f v_max, float R);
	virtual ~VelocityGrid();
	const DeviceVelocityGrid& device() const;
};

class DeviceVelocityGrid: public VelocityGridBase
{
	friend class VelocityGrid;
	DeviceVelocityGrid(const VelocityGrid&);
public:
//	void free();
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

