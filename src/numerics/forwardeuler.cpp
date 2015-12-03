/*
 * forwardeuler.cpp
 *
 *  Created on: Nov 13, 2014
 *      Author: tony
 */

#include "forwardeuler.h"

namespace solvers {

bool forward_euler::step( const double told, const double tnew, double* Yold, double* Ynew ){
	cells->celldiffeq( told, Yold, dYdt);

	// Compute the new states for each cell
	for (int i=0; i< neqns; i++)
		Ynew[i] = Yold[i] + dYdt[i] * (tnew - told); // Forward Euler

	return true; // successful step

}

} /* namespace solvers */
