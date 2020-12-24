/*
 * collision_integral.h
 *
 *  Created on: 17.05.2010
 *      Author: serge
 */

#include "options.hpp"

#ifndef COLLISION_INTEGRAL_H_
void evolve_colisions (float * f,
		       float * direct_integral,
		       float * inverse_integral,
		       float * b,
		       float * a,
		       float * correction_array,
			   const Options& opts);

void relax_f (float * f,
	      float * frequency,
	      float * inverse_integral,
	      float courant,
	      float * time,
		  const Options& opts);
#define COLLISION_INTEGRAL_H_
#endif /* COLLISION_INTEGRAL_H_ */
