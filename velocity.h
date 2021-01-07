/*
 * velocity.h
 *
 *  Created on: 17.02.2010
 *      Author: serge
 */

#ifndef VELOCITY_H_
#define VELOCITY_H_

#include "common.hpp"
#include "velocity_grid.hpp"
#include "options.hpp"

#define N_X 16
#define N_Y 16
#define N_Z 16
#define N_YZ (N_Y*N_Z) 

#define NP_X 256

#define XI_X(v)    (v)
#define XI_Y(v)    (&((v)[opts.nxyz]))
#define XI_Z(v)    (&((v)[2 * opts.nxyz]))
#define XI_2(v)    (&((v)[3 * opts.nxyz]))
#define XI_X_X(v)  (&((v)[4 * opts.nxyz]))
#define XI_X_Y(v)  (&((v)[5 * opts.nxyz]))
#define XI_X_Z(v)  (&((v)[6 * opts.nxyz]))
#define XI_X_2(v)  (&((v)[7 * opts.nxyz]))
#define XI_Y_Y(v)  (&((v)[8 * opts.nxyz]))
#define XI_Y_Z(v)  (&((v)[9 * opts.nxyz]))
#define XI_Y_2(v)  (&((v)[10 * opts.nxyz]))
#define XI_Z_Z(v)  (&((v)[11 * opts.nxyz]))
#define XI_Z_2(v)  (&((v)[12 * opts.nxyz]))
#define XI_2_2(v)  (&((v)[13 * opts.nxyz]))

#define ERR 1.0e-10

void free_g_params (int ** h_g_params, int ** d_g_params);

void init_correction_array (float ** correction_array, const DeviceVelocityGrid&, const Options& opts);

#endif /* VELOCITY_H_ */
