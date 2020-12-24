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

typedef struct _IntVector IntVector;
typedef struct _DVector DVector;

struct _IntVector {
  int x, y, z;
};

struct _DVector {
  float x, y, z;
};

#endif // COMMON_HPP
