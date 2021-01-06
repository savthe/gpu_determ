#ifndef COMMON_HPP
#define COMMON_HPP

struct IntParams
{
	int n_pnt;
	int n_u;
	int n_v;
	int n_w;
	int n_uvw;
};

/*
typedef struct _IntVector IntVector;
typedef struct _DVector DVector;

struct _IntVector {
  int x, y, z;
};

struct _DVector {
  float x, y, z;
};
*/

template<typename T> 
struct Vector3D
{
	Vector3D(): x(), y(), z() {}
	Vector3D(const T& x_, const T& y_, const T& z_): x(x_), y(y_), z(z_) {}
	T x, y, z;
};

using Vector3i = Vector3D<int>;
using Vector3f = Vector3D<float>;
#endif // COMMON_HPP
