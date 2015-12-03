/*
 * rk45.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: tony
 */

#include "rk45.h"
#include <cmath>

namespace solvers {

bool rk45::step( const double told, const double tnew, double* Yold, double* Ynew ){
	const double dt  = tnew - told;

	// First function call:		k1 = dt * f(tn, yn)
	cells->celldiffeq( told, Yold, k1);
	// Finish computing k1 and pre-compute next y
	for (int i=0; i< neqns ; i++) {
		k1[i] = k1[i] * dt;
		Yhalf[i] = Yold[i] + c21 * k1[i];
	}

	// Second function call:	k2 = dt * f(t + c20 * dt, yn + c21 * k1)
	const double tpc20 = told + c20*dt;
	cells->celldiffeq( tpc20, Yhalf, k2);
	for (int i=0; i< neqns ; i++){	// Finish computing k2 and pre-compute next y
		k2[i] = k2[i] * dt;
		Yhalf[i] = Yold[i] + c31 * k1[i] + c32 * k2[i];
	}

	// Third function call		k3 = dt * f(t + c30 * dt, yn + c31 * k1 + c32 * k2)
	const double tpc30 = told + c30*dt;
	cells->celldiffeq( tpc30, Yhalf, k3);
	for (int i=0; i< neqns ; i++){	// Finish computing k3 and pre-compute next y
		k3[i] = k3[i] * dt;
		Yhalf[i] = Yold[i] + c41 * k1[i] + c42 * k2[i] + c43 * k3[i];
	}

	// Fourth function call		k4 = dt * f(t + c40 * dt, yn + c41 * k1 + c42 * k2 + c43 * k3)
	const double tpc40 = told + c40*dt;
	cells->celldiffeq( tpc40, Yhalf, k4);
	for (int i=0; i< neqns ; i++){	// Finish computing k4 and pre-compute next y
		k4[i] = k4[i] * dt;
		Yhalf[i] = Yold[i] + c51 * k1[i] + c52 * k2[i] + c53 * k3[i] + c54 * k4[i];
	}

	// Fifth function call		k5 = dt * f(t + dt, yn + c51 * k1 + c52 * k2 + c53 * k3 + c54 * k4 )
	const double tpc50 = tnew;
	cells->celldiffeq( tpc50, Yhalf, k5);
	for (int i=0; i< neqns ; i++){	// Finish computing k5 and pre-compute next y and the complete 4th order approximation
		k5[i] = k5[i] * dt;
		Yhalf[i] = Yold[i] + c61 * k1[i] + c62 * k2[i] + c63 * k3[i] + c64 * k4[i] + c65 * k5[i];
		//  y4 = yn + a1 * k1 + a3 * k3 + a4 * k4 + a5 * k5;
		Y4[i]   = Yold[i] + a1 * k1[i] + a3 * k3[i] + a4 * k4[i] + a5 * k5[i];
	}

	// Sixth function call		k6 = dt * f(t + c60 * dt, yn + c61 * k1 + c62 * k2 + c63 * k3 + c64 * k4 + c65 * k5);
	const double tpc60 = told + c60*dt;
	cells->celldiffeq( tpc60, Yhalf, k6);
	epsilon = 0; // diff between 4th and 5th order results
	for (int i=0; i< neqns ; i++){	// Finish computing k5, and the complete 5th order approximation
		k6[i] = k6[i] * dt;
		//  y5 = yn + b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6;
		Ynew[i] = Yold[i] + b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i];
		//  epsilon = fabs(y4 - y5);
		epsilon = epsilon + fabs(Ynew[i] - Y4[i]);
	}
	epsilon = epsilon/( neqns );

	return true; // successful step
}

} /* namespace solvers */
