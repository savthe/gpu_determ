/*
 * compute_matrix.cu
 *
 *  Created on: 17.02.2010
 *      Author: serge
 */

//#include <math.h>
#define _USE_MATH_DEFINES
#include <assert.h>
#include <cuda.h>
#include <stdio.h>
#include "velocity.h"
#include "compute_matrix.h"
#include "gputimer.hpp"
#include <math.h>


#if 1

__device__ double f1 (double p, double A)
{
  
  double sqap;

  sqap = 1. - A*A - p*p;
  
  if ((A < ERR) && (A > -ERR))
    return M_PI_2 * p;
  else if ((p < ERR) && (p > -ERR))
    return M_PI_2 * A;
  else if ((sqap < ERR) && (sqap > -ERR))
    return M_PI_2;
  else if ((sqap >= ERR) && (p >= ERR) && (A >= ERR)) {
    sqap = sqrt (sqap);
    return p * atan (sqap/A) + A * atan (sqap/p) + atan (A*p/sqap);
  }
  //  else
    //    assert (0);

  return 0.;
}


__device__ double f2 (double p, double A)
{
  double sqap;
  
  sqap = 1. - A*A -p*p;
  
  if (A > (1 - ERR))
    return 0.;
  else if ((p < ERR) && (p > -ERR))
    return M_PI_2 * A;
  else if ((sqap < ERR) && (sqap > -ERR))
    return M_PI_2;
  else if ((sqap >= ERR) && (p >= ERR)) {
    sqap = sqrt (sqap);
    return p * atan (sqap/A) + A * atan (sqap/p) + atan (A*p/sqap);
  }
  //  else
  //    assert (0);
  return 0.;
}


__device__ double f3 (double p, double A)
{
  double sqap;
  
  sqap = 1. - A*A - p*p;
  
  if (A < ERR)
    return 0.;
  else if ((p < ERR) && (p > -ERR))
    return -M_PI_2 * A;
  else if ((sqap < ERR) && (sqap > -ERR))
    return M_PI_2 * (p - 1.);
  else if ((sqap >= ERR) && (p >= ERR)) {
    sqap = sqrt (sqap);
    return M_PI_2 * p - p * atan (sqap/A) - A * atan (sqap/p) - atan (A*p/sqap);
  }
  //  else
  //    assert (0);
  return 0.;
}


__device__ double f4 (double p, double A)
{
  double sqap;
  
  sqap = 1. - A*A - p*p;
  
  if (A > (1.-ERR))
    return M_PI_2 * p;
  else if ((p < ERR) && (p > -ERR))
    return -M_PI_2 * A;
  else if((sqap < ERR) && (sqap > -ERR))
    return M_PI_2 * (p - 1.);
  else if((sqap >= ERR) && (p >= ERR)) {
    sqap = sqrt (sqap);
    return M_PI_2 * p - p * atan (sqap/A) - A * atan (sqap/p) - atan (A*p/sqap);
  }
  //  else
  //   assert (0);
  return 0.;
}

__device__ double calc_matrix_element (double du, double dv, double dw, 
				       int in, int jn, int kn, 
				       double g)
{
  double u1, u2, v1, v2, w1, w2, sqg = g*g;
  double uu1, uu2, vv1, vv2, ww1, ww2;
  double eps1, eps3;
  double sqp[4], p[4], ppz[2];
  int i;
  double A1, A2, A3, A4;

  u1 = in > 0 ? du * (in - 1) : 0.;
  u2 = du * (in + 1);
  v1 = jn > 0 ? dv * (jn - 1) : 0.;
  v2 = dv * (jn + 1);
  w1 = kn > 0 ? dw * (kn - 1) : 0.;
  w2 = dw * (kn + 1);
  
  uu1 = u1 * u1;
  uu2 = u2 * u2;
  vv1 = v1 * v1;
  vv2 = v2 * v2;
  ww1 = w1 * w1;
  ww2 = w2 * w2;
  
  if (uu1 + vv1 + ww1 < sqg && uu2 + vv2 + ww2 > sqg) {
    if (in > 1)
      eps1 = atan (sqrt (sqg/uu1 - 1.));
    else
      eps1 = M_PI_2;
    
    if (jn > 1)
      eps3 = M_PI_2 - atan (sqrt (sqg/vv1 - 1.));
    else
      eps3 = 0;
  }
  else {
    eps1 = 0;
    eps3 = 1;
  }
  
  if (eps1 > eps3) {
    sqp[0] = 1. - (uu2 + vv2)/sqg;
    sqp[1] = 1. - (uu1 + vv2)/sqg;
    sqp[2] = 1. - (uu1 + vv1)/sqg;
    sqp[3] = 1. - (uu2 + vv1)/sqg;
    
    for (i = 0; i < 4; i++)
      p[i] = sqp[i] < 0 ? 0. : sqrt (sqp[i]);
    ppz[0] = w1/g;
    ppz[1] = w2/g > 1. ? 1. : w2/g;
    
    A1 = u1/g;
    A2 = u2/g;
    A3 = v1/g;
    A4 = v2/g;
    
    if (p[1] < p[3])
      if (ppz[0] <= p[0])
	if (ppz[1] <= p[0])
	  return 0.;
	else if (ppz[1] <= p[1])
	  return f4(ppz[1], A4) - f4(p[0], A4) - f2(ppz[1], A2) + f2(p[0], A2);
	else if (ppz[1] <= p[3])
	  return f4(p[1], A4) - f4(p[0], A4) + f1 (ppz[1], A1) - f1(p[1], A1) 
	    - f2(ppz[1], A2) + f2(p[0], A2);
	else if (ppz[1] <= p[2])
	  return f4(p[1], A4) - f4(p[0], A4) + f1(ppz[1], A1) - f1(p[1], A1)-
	    f2(p[3], A2) + f2(p[0], A2) - f3(ppz[1], A3) + f3(p[3], A3);
	else
	  return f4(p[1], A4) - f4(p[0], A4) + f1(p[2], A1) - f1(p[1], A1) -
	    f2(p[3], A2) + f2(p[0], A2) - f3(p[2], A3) + f3(p[3], A3);
      else if (ppz[0] <= p[1])
	if (ppz[1] <= p[1])
	  return f4(ppz[1], A4) - f4(ppz[0], A4) - f2(ppz[1], A2) + f2(ppz[0], A2);
	else if (ppz[1] <= p[3])
	  return f4(p[1], A4) - f4(ppz[0], A4) + f1(ppz[1], A1) - f1(p[1], A1) -
	    f2(ppz[1], A2) + f2(ppz[0], A2);
	else if (ppz[1] <= p[2])
	  return f4(p[1], A4) - f4(ppz[0], A4) + f1(ppz[1], A1) - f1(p[1], A1) -
	    f2(p[3], A2) + f2(ppz[0], A2) - f3(ppz[1], A3) + f3(p[3], A3);
	else
	  return f4(p[1], A4) - f4(ppz[0], A4) + f1(p[2], A1) - f1(p[1], A1) -
	    f2(p[3], A2) + f2(ppz[0], A2) - f3(p[2], A3) + f3(p[3], A3);
      else if (ppz[0]<=p[3])
	if (ppz[1]<=p[3])
	  return f1(ppz[1], A1) - f1(ppz[0], A1) - f2(ppz[1], A2) + f2(ppz[0], A2);
	else if (ppz[1]<=p[2])
	  return f1(ppz[1], A1) - f1(ppz[0], A1) - f2(p[3], A2) + f2(ppz[0], A2) - 
	    f3(ppz[1], A3) + f3(p[3], A3);
	else
	  return f1(p[2], A1) - f1(ppz[0], A1) - f2(p[3], A2) + f2(ppz[0], A2) - 
	    f3(p[2], A3) + f3(p[3], A3);
      else if (ppz[0]<=p[2])
	if (ppz[1]<=p[2])
	  return f1(ppz[1], A1) - f1(ppz[0], A1) - f3(ppz[1], A3) + f3(ppz[0], A3);
	else
	  return f1(p[2], A1) - f1(ppz[0], A1) - f3(p[2], A3) + f3(ppz[0], A3);
      else
	return 0.;
    else {
      if (ppz[0]<=p[0])
	if (ppz[1]<=p[0])
	  return 0.;
	else if (ppz[1]<=p[3])
	  return f4(ppz[1], A4) - f4(p[0], A4) - f2(ppz[1], A2) + f2(p[0], A2);
	else if (ppz[1]<=p[1])
	  return f4(ppz[1], A4) - f4(p[0], A4) - f2(p[3], A2) + f2(p[0], A2) - 
	    f3(ppz[1], A3) + f3(p[3], A3);
	else if (ppz[1]<=p[2])
	  return f4(p[1], A4) - f4(p[0], A4) + f1(ppz[1], A1) - f1(p[1], A1) -
	    f2(p[3], A2) + f2(p[0], A2) - f3(ppz[1], A3) + f3(p[3], A3);
	else
	  return f4(p[1], A4) - f4(p[0], A4) + f1(p[2], A1) - f1(p[1], A1) -
	    f2(p[3], A2) + f2(p[0], A2) - f3(p[2], A3) + f3(p[3], A3);
      else if (ppz[0]<=p[3])
	if (ppz[1]<=p[3])
	  return f4(ppz[1], A4) - f4(ppz[0], A4) - f2(ppz[1], A2) + f2(ppz[0], A2);
	else if (ppz[1]<=p[1])
	  return f4(ppz[1], A4) - f4(ppz[0], A4) - f2(p[3], A2) + f2(ppz[0], A2) - 
	    f3(ppz[1], A3) + f3(p[3], A3);
	else if (ppz[1]<=p[2])
	  return f4(p[1], A4) - f4(ppz[0], A4) + f1(ppz[1], A1) - f1(p[1], A1) -
	    f2(p[3], A2) + f2(ppz[0], A2) - f3(ppz[1], A3) + f3(p[3], A3);
	else
	  return f4(p[1], A4) - f4(ppz[0], A4) + f1(p[2], A1) - f1(p[1], A1) -
	    f2(p[3], A2) + f2(ppz[0], A2) - f3(p[2], A3) + f3(p[3], A3);
      else if (ppz[0]<=p[1])
	if (ppz[1]<=p[1])
	  return f4(ppz[1], A4) - f4(ppz[0], A4) - f3(ppz[1], A3) + f3(ppz[0], A3);
	else if (ppz[1]<=p[2])
	  return f4(p[1], A4) - f4(ppz[0], A4) + f1(ppz[1], A1) - f1(p[1], A1) - 
	    f3(ppz[1], A3) + f3(ppz[0], A3);
	else
	  return f4(p[1], A4) - f4(ppz[0], A4) + f1(p[2], A1) - f1(p[1], A1) - 
	    f3(p[2], A3) + f3(ppz[0], A3);
      else if (ppz[0]<=p[2])
	if (ppz[1]<=p[2])
	  return f1(ppz[1], A1) - f1(ppz[0], A1) - f3(ppz[1], A3) + f3(ppz[0], A3);
	else
	  return f1(p[2], A1) - f1(ppz[0], A1) - f3(p[2], A3) + f3(ppz[0], A3);
      else
	return 0.;
    }
  }
  else
    return 0.;
}

#else

__device__ float f1 (float p, float A)
{
  
  float sqap;

  sqap = 1. - A*A - p*p;
  
  if ((A < ERR) && (A > -ERR))
    return M_PI_2 * p;
  else if ((p < ERR) && (p > -ERR))
    return M_PI_2 * A;
  else if ((sqap < ERR) && (sqap > -ERR))
    return M_PI_2;
  else if ((sqap >= ERR) && (p >= ERR) && (A >= ERR)) {
    sqap = sqrt (sqap);
    return p * atan (sqap/A) + A * atan (sqap/p) + atan (A*p/sqap);
  }
  //  else
    //    assert (0);

  return 0.;
}


__device__ float f2 (float p, float A)
{
  float sqap;
  
  sqap = 1. - A*A -p*p;
  
  if (A > (1 - ERR))
    return 0.;
  else if ((p < ERR) && (p > -ERR))
    return M_PI_2 * A;
  else if ((sqap < ERR) && (sqap > -ERR))
    return M_PI_2;
  else if ((sqap >= ERR) && (p >= ERR)) {
    sqap = sqrt (sqap);
    return p * atan (sqap/A) + A * atan (sqap/p) + atan (A*p/sqap);
  }
  //  else
  //    assert (0);
  return 0.;
}


__device__ float f3 (float p, float A)
{
  float sqap;
  
  sqap = 1. - A*A - p*p;
  
  if (A < ERR)
    return 0.;
  else if ((p < ERR) && (p > -ERR))
    return -M_PI_2 * A;
  else if ((sqap < ERR) && (sqap > -ERR))
    return M_PI_2 * (p - 1.);
  else if ((sqap >= ERR) && (p >= ERR)) {
    sqap = sqrt (sqap);
    return M_PI_2 * p - p * atan (sqap/A) - A * atan (sqap/p) - atan (A*p/sqap);
  }
  //  else
  //    assert (0);
  return 0.;
}


__device__ float f4 (float p, float A)
{
  float sqap;
  
  sqap = 1. - A*A - p*p;
  
  if (A > (1.-ERR))
    return M_PI_2 * p;
  else if ((p < ERR) && (p > -ERR))
    return -M_PI_2 * A;
  else if((sqap < ERR) && (sqap > -ERR))
    return M_PI_2 * (p - 1.);
  else if((sqap >= ERR) && (p >= ERR)) {
    sqap = sqrt (sqap);
    return M_PI_2 * p - p * atan (sqap/A) - A * atan (sqap/p) - atan (A*p/sqap);
  }
  //  else
  //   assert (0);
  return 0.;
}

__device__ float calc_matrix_element (float * float_params, 
				       int in, int jn, int kn, 
				       float g)
{
  float u1, u2, v1, v2, w1, w2, sqg = g*g;
  float uu1, uu2, vv1, vv2, ww1, ww2;
  float eps1, eps3;
  float sqp[4], p[4], ppz[2];
  int i;
  float A1, A2, A3, A4;
  float du, dv, dw;

  du  = D_U (float_params);
  dv  = D_V (float_params);
  dw  = D_W (float_params);

  u1 = in > 0 ? du * (in - 1) : 0.;
  u2 = du * (in + 1);
  v1 = jn > 0 ? dv * (jn - 1) : 0.;
  v2 = dv * (jn + 1);
  w1 = kn > 0 ? dw * (kn - 1) : 0.;
  w2 = dw * (kn + 1);
  
  uu1 = u1 * u1;
  uu2 = u2 * u2;
  vv1 = v1 * v1;
  vv2 = v2 * v2;
  ww1 = w1 * w1;
  ww2 = w2 * w2;
  
  if (uu1 + vv1 + ww1 < sqg && uu2 + vv2 + ww2 > sqg) {
    if (in > 1)
      eps1 = atan (sqrt (sqg/uu1 - 1.));
    else
      eps1 = M_PI_2;
    
    if (jn > 1)
      eps3 = M_PI_2 - atan (sqrt (sqg/vv1 - 1.));
    else
      eps3 = 0;
  }
  else {
    eps1 = 0;
    eps3 = 1;
  }
  
  if (eps1 > eps3) {
    sqp[0] = 1. - (uu2 + vv2)/sqg;
    sqp[1] = 1. - (uu1 + vv2)/sqg;
    sqp[2] = 1. - (uu1 + vv1)/sqg;
    sqp[3] = 1. - (uu2 + vv1)/sqg;
    
    for (i = 0; i < 4; i++)
      p[i] = sqp[i] < 0 ? 0. : sqrt (sqp[i]);
    ppz[0] = w1/g;
    ppz[1] = w2/g > 1. ? 1. : w2/g;
    
    A1 = u1/g;
    A2 = u2/g;
    A3 = v1/g;
    A4 = v2/g;
    
    if (p[1] < p[3])
      if (ppz[0] <= p[0])
	if (ppz[1] <= p[0])
	  return 0.;
	else if (ppz[1] <= p[1])
	  return f4(ppz[1], A4) - f4(p[0], A4) - f2(ppz[1], A2) + f2(p[0], A2);
	else if (ppz[1] <= p[3])
	  return f4(p[1], A4) - f4(p[0], A4) + f1 (ppz[1], A1) - f1(p[1], A1) 
	    - f2(ppz[1], A2) + f2(p[0], A2);
	else if (ppz[1] <= p[2])
	  return f4(p[1], A4) - f4(p[0], A4) + f1(ppz[1], A1) - f1(p[1], A1)-
	    f2(p[3], A2) + f2(p[0], A2) - f3(ppz[1], A3) + f3(p[3], A3);
	else
	  return f4(p[1], A4) - f4(p[0], A4) + f1(p[2], A1) - f1(p[1], A1) -
	    f2(p[3], A2) + f2(p[0], A2) - f3(p[2], A3) + f3(p[3], A3);
      else if (ppz[0] <= p[1])
	if (ppz[1] <= p[1])
	  return f4(ppz[1], A4) - f4(ppz[0], A4) - f2(ppz[1], A2) + f2(ppz[0], A2);
	else if (ppz[1] <= p[3])
	  return f4(p[1], A4) - f4(ppz[0], A4) + f1(ppz[1], A1) - f1(p[1], A1) -
	    f2(ppz[1], A2) + f2(ppz[0], A2);
	else if (ppz[1] <= p[2])
	  return f4(p[1], A4) - f4(ppz[0], A4) + f1(ppz[1], A1) - f1(p[1], A1) -
	    f2(p[3], A2) + f2(ppz[0], A2) - f3(ppz[1], A3) + f3(p[3], A3);
	else
	  return f4(p[1], A4) - f4(ppz[0], A4) + f1(p[2], A1) - f1(p[1], A1) -
	    f2(p[3], A2) + f2(ppz[0], A2) - f3(p[2], A3) + f3(p[3], A3);
      else if (ppz[0]<=p[3])
	if (ppz[1]<=p[3])
	  return f1(ppz[1], A1) - f1(ppz[0], A1) - f2(ppz[1], A2) + f2(ppz[0], A2);
	else if (ppz[1]<=p[2])
	  return f1(ppz[1], A1) - f1(ppz[0], A1) - f2(p[3], A2) + f2(ppz[0], A2) - 
	    f3(ppz[1], A3) + f3(p[3], A3);
	else
	  return f1(p[2], A1) - f1(ppz[0], A1) - f2(p[3], A2) + f2(ppz[0], A2) - 
	    f3(p[2], A3) + f3(p[3], A3);
      else if (ppz[0]<=p[2])
	if (ppz[1]<=p[2])
	  return f1(ppz[1], A1) - f1(ppz[0], A1) - f3(ppz[1], A3) + f3(ppz[0], A3);
	else
	  return f1(p[2], A1) - f1(ppz[0], A1) - f3(p[2], A3) + f3(ppz[0], A3);
      else
	return 0.;
    else {
      if (ppz[0]<=p[0])
	if (ppz[1]<=p[0])
	  return 0.;
	else if (ppz[1]<=p[3])
	  return f4(ppz[1], A4) - f4(p[0], A4) - f2(ppz[1], A2) + f2(p[0], A2);
	else if (ppz[1]<=p[1])
	  return f4(ppz[1], A4) - f4(p[0], A4) - f2(p[3], A2) + f2(p[0], A2) - 
	    f3(ppz[1], A3) + f3(p[3], A3);
	else if (ppz[1]<=p[2])
	  return f4(p[1], A4) - f4(p[0], A4) + f1(ppz[1], A1) - f1(p[1], A1) -
	    f2(p[3], A2) + f2(p[0], A2) - f3(ppz[1], A3) + f3(p[3], A3);
	else
	  return f4(p[1], A4) - f4(p[0], A4) + f1(p[2], A1) - f1(p[1], A1) -
	    f2(p[3], A2) + f2(p[0], A2) - f3(p[2], A3) + f3(p[3], A3);
      else if (ppz[0]<=p[3])
	if (ppz[1]<=p[3])
	  return f4(ppz[1], A4) - f4(ppz[0], A4) - f2(ppz[1], A2) + f2(ppz[0], A2);
	else if (ppz[1]<=p[1])
	  return f4(ppz[1], A4) - f4(ppz[0], A4) - f2(p[3], A2) + f2(ppz[0], A2) - 
	    f3(ppz[1], A3) + f3(p[3], A3);
	else if (ppz[1]<=p[2])
	  return f4(p[1], A4) - f4(ppz[0], A4) + f1(ppz[1], A1) - f1(p[1], A1) -
	    f2(p[3], A2) + f2(ppz[0], A2) - f3(ppz[1], A3) + f3(p[3], A3);
	else
	  return f4(p[1], A4) - f4(ppz[0], A4) + f1(p[2], A1) - f1(p[1], A1) -
	    f2(p[3], A2) + f2(ppz[0], A2) - f3(p[2], A3) + f3(p[3], A3);
      else if (ppz[0]<=p[1])
	if (ppz[1]<=p[1])
	  return f4(ppz[1], A4) - f4(ppz[0], A4) - f3(ppz[1], A3) + f3(ppz[0], A3);
	else if (ppz[1]<=p[2])
	  return f4(p[1], A4) - f4(ppz[0], A4) + f1(ppz[1], A1) - f1(p[1], A1) - 
	    f3(ppz[1], A3) + f3(ppz[0], A3);
	else
	  return f4(p[1], A4) - f4(ppz[0], A4) + f1(p[2], A1) - f1(p[1], A1) - 
	    f3(p[2], A3) + f3(ppz[0], A3);
      else if (ppz[0]<=p[2])
	if (ppz[1]<=p[2])
	  return f1(ppz[1], A1) - f1(ppz[0], A1) - f3(ppz[1], A3) + f3(ppz[0], A3);
	else
	  return f1(p[2], A1) - f1(ppz[0], A1) - f3(p[2], A3) + f3(ppz[0], A3);
      else
	return 0.;
    }
  }
  else
    return 0.;
}
#endif

__global__ void store_matrix_elements (const VelocityGrid vgrid, float * b, float * a)
{
  short int index, index1;

  //float* float_params = vgrid.float_params;
  
  float g, gu, gv, gw;
  short int gi, gj, gk, in, jn, kn, jn1, kn1;
  
  if (blockIdx.x >= N_X) return;

  gi = blockIdx.x;
  gu = vgrid.d_u * gi;

  if (blockIdx.y == 0) 
    for (index = threadIdx.x; index < N_YZ; index += blockDim.x) {
      gj = index / N_Z;
      gk = index % N_Z;

      gv = vgrid.d_v * gj;
      gw = vgrid.d_w * gk;
      
      g = sqrt (gu*gu + gv*gv + gw*gw);

      b[gi * N_YZ + index] = 4.0 * M_PI * g * vgrid.d3v;
    }

  if (blockIdx.y >= (N_X >> 1)) return;

  in = (gi & 1) | (blockIdx.y << 1);

  for (jn1 = 0, index = 0; jn1 < (N_Y >> 1); jn1++) 
    for (kn1 = 0; kn1 < (N_Z >> 1); kn1++, index++) { 
      for (index1 = threadIdx.x; index1 < N_YZ; index1 += blockDim.x) {
	gj = index1 / N_Z;
	gk = index1 % N_Z;

	double du = vgrid.d_u;
	double dv = vgrid.d_v;
	double dw = vgrid.d_w;

	jn = (gj & 1) | (jn1 << 1);
	kn = (gk & 1) | (kn1 << 1);

	gv = vgrid.d_v * gj;
	gw = vgrid.d_w * gk;
	
	//	g = sqrt (gu*gu + gv*gv + gw*gw);
	double g1 = sqrt (gu*gu + gv*gv + gw*gw);

	double A = calc_matrix_element (du, dv, dw, in, jn, kn, g1);
	a[(((in>>1) * N_X + gi) * (N_YZ >> 2) + 
	   index) * N_YZ + index1] = 
	  (A < ERR) ? 0.0f : A * 2.0f * g1 * vgrid.d3v;
	//	    (A < ERR) ? 0.0f : A;
      } 
    }
}

void init_matrices (const VelocityGrid& vgrid, float ** b, float ** a, const Options& opts)
{
	
  	GpuTimer timer;
  	timer.start();
  dim3 dg, db;
  
  dg.x = N_X;
  dg.y = N_X >> 1;
  dg.z = 1;
  db.x = 128;
  db.y = 1;
  db.z = 1;


  cudaMalloc ((void **) a, 
	      opts.nxyz * (opts.nxyz >> 3) * sizeof (float));

  printf ("ALLOCATED %d BYTES, FOR INVERSE COLISIONS MATRIX: %s\n", 
	  opts.nxyz * (opts.nxyz >> 3) *  sizeof (float),
	  cudaGetErrorString (cudaGetLastError ()));

  cudaMalloc ((void **) b, opts.nxyz * sizeof (float));

  printf ("ALLOCATED %d BYTES, FOR DIRECT COLISIONS: %s\n", 
	  opts.nxyz * sizeof (float),
	  cudaGetErrorString (cudaGetLastError ()));

  store_matrix_elements <<<dg, db>>> (vgrid, *b, *a); 
  printf ("MATRICES STORED IN GLOBAL MEMORY: %s\n", 
	  cudaGetErrorString (cudaGetLastError ()));
  
  cudaThreadSynchronize();

  	timer.stop();
	printf ("MATRICES INITIALIZATION TOOK %f MS!\n", timer.elapsed());
}
