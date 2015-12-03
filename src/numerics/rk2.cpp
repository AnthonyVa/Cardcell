/*
 * rk2.cpp
 *
 *  Created on: Nov 13, 2014
 *      Author: tony
 */

#include "rk2.h"

namespace solvers {


bool rk2::step( const double told, const double tnew, double* Yold, double* Ynew ){
	const double dt  = tnew - told;
	const double dt2 = dt/2;
	const double thalf = told + dt2;

	// First function call			k1 = dt f(tn, yn) / 2
	cells->celldiffeq( told, Yold, k1);
	// Finish computing k1 and pre-compute Yhalf for next step
	for (int i=0; i< neqns; i++) {
		k1[i] = k1[i] * dt2;
		// Yhalf is (yn + k1)
		Yhalf[i] = Yold[i] + k1[i];
	}

	// Second function call			k2 = dt f(tn+1/2, yn + k1)
	cells->celldiffeq( thalf, Yhalf, dYdt);
	// Compute the new state
	for (int i=0; i< neqns; i++)
		Ynew[i] = Yold[i] + dYdt[i] * dt; //

	return true; // successful step
}

} /* namespace solvers */
