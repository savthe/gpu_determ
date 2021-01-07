#include <cassert>
#include "velocity_grid.hpp"
#include "velocity.h"


bool check_if_point_is_in_domain (Vector3i n_points, Vector3f v_min, Vector3f v_max,
			                     int i, int j, int k, float R)
{
	float u = (v_max.x - v_min.x) * (2*i + 1 - n_points.x) / 2. / n_points.x; 
	float v = (v_max.y - v_min.y) * (2*j + 1 - n_points.y) / 2. / n_points.y;
	float w = (v_max.z - v_min.z) * (2*k + 1 - n_points.z) / 2. / n_points.z; 

	return R <= 0 || u*u + v*v + w*w < R*R;
}


VelocityGrid::VelocityGrid(Vector3i n_points, Vector3f v_min, Vector3f v_max, float r)
{
	init(n_points, v_min, v_max, r);
}

const DeviceVelocityGrid& VelocityGrid::device() const
{
	static DeviceVelocityGrid g(*this);
	return g;
}

DeviceVelocityGrid::DeviceVelocityGrid(const VelocityGrid& g):
	VelocityGridBase(g)
{
	cudaMalloc ((void **) &u_index, n_pnt * sizeof (float));
	cudaMemcpy (u_index, g.u_index, n_pnt * sizeof(float), cudaMemcpyHostToDevice); 

	cudaMalloc ((void **) &v_index, n_pnt * sizeof (float));
	cudaMemcpy (v_index, g.v_index, n_pnt * sizeof(float), cudaMemcpyHostToDevice); 

	cudaMalloc ((void **) &w_index, n_pnt * sizeof (float));
	cudaMemcpy (w_index, g.w_index, n_pnt * sizeof(float), cudaMemcpyHostToDevice); 

	cudaMalloc ((void **) &u, n_u * sizeof (float));
	cudaMemcpy (u, g.u, n_u * sizeof(float), cudaMemcpyHostToDevice); 

	cudaMalloc ((void **) &w, n_w * sizeof (float));
	cudaMemcpy (w, g.w, n_w * sizeof(float), cudaMemcpyHostToDevice); 

	cudaMalloc ((void **) &v, n_v * sizeof (float));
	cudaMemcpy (v, g.v, n_v * sizeof(float), cudaMemcpyHostToDevice); 

}

VelocityGrid::~VelocityGrid()
{
	if(!u) return;
	delete[] u;
	delete[] v;
	delete[] w;
	delete[] u_index;
	delete[] v_index;
	delete[] w_index;
}

void VelocityGrid::init(Vector3i n_points, Vector3f v_min, Vector3f v_max, float R)
{
  
	assert (n_points.x > 0 && n_points.y > 0 && n_points.z > 0);
	assert (v_min.x < v_max.x && v_min.y < v_max.y && v_min.z < v_max.z);
 
	n_pnt = 0;
	for (int i = 0; i < n_points.x; i++)
	for (int j = 0; j < n_points.y; j++)
	for (int k = 0; k < n_points.z; k++)
		n_pnt += check_if_point_is_in_domain (n_points, v_min, v_max, i, j, k, R);
  
	n_u = n_points.x;
	n_v = n_points.y;
	n_w = n_points.z;
	n_uvw = n_points.x * n_points.y * n_points.z;

	u = new float[n_u];
	v = new float[n_v];
	w = new float[n_w];

	u_index = new int[n_pnt];
	v_index = new int[n_pnt];
	w_index = new int[n_pnt];

	min_u = v_min.x;
	min_v = v_min.y;
	min_w = v_min.z;
	max_u = v_max.x;
	max_v = v_max.y;
	max_w = v_max.z;

	d_u = (max_u - min_u) / n_u;
	d_v = (max_v - min_v) / n_v;
	d_w = (max_w - min_w) / n_w;

	d3v = d_u*d_v*d_w;
  
  
	for (int i = 0; i < n_u; i++) u[i] = min_u + d_u * (2.*i+1.)/2.;
	for (int i = 0; i < n_v; i++) v[i] = min_v + d_v * (2.*i+1.)/2.;
	for (int i = 0; i < n_w; i++) w[i] = min_w + d_w * (2.*i+1.)/2.;

	int index = 0;
	for (int i = 0; i < n_u; i++)
	for (int j = 0; j < n_v; j++)
	for (int k = 0; k < n_w; k++)
		if (check_if_point_is_in_domain (n_points, v_min, v_max, i, j, k, R)) 
		{
			u_index[index] = i;
			v_index[index] = j;
			w_index[index] = k;
			index++;
	    }
}

