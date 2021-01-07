/*
 * compute_matrix.h
 *
 *  Created on: 17.05.2010
 *      Author: serge
 */

#ifndef COMPUTE_MATRIX_H_

#include "options.hpp"
#include "velocity_grid.hpp"

void init_matrices (const DeviceVelocityGrid&, float ** b, float ** a, const Options& opts);
#define COMPUTE_MATRIX_H_
#endif /* COMPUTE_MATRIX_H_ */
