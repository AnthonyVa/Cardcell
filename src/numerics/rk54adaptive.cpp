/*
 * rk54adaptive.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: W1084063
 */

#include "rk54adaptive.h"
#include <algorithm>    // std::max

namespace solvers {

/*
 * Adaptive time stepping:
 * 1. Start at time told and starting value Yold; set the step size to be dt
 * 2. Use rk45 to compute the approximation Ynew and error estimate epsilon
 * 3. check error and adjust step size
 */
bool rk54adaptive::step( const double told, const double tnew, double* Yold, double* Ynew ){
	double t = told;
	double step = (tnew - told)/2;
	int integrationFlag = 1; // 0 == good, 1 == not good or time to quit
	int badStep = 0; // count the number of bad steps and quit if too many

	for (int k=0; k < maxIterations; k++) {
		if (fabs(step ) < minStep)
			step = sgn(step) * minStep;
		if (fabs(step) > maxStep)
			step = sgn(step) * maxStep;

		const double d = fabs(tnew - t);
		if (d <= fabs(step) ) {
			integrationFlag = 0; // 0 == good, 1 == time to quit
			if (d <= ( delta * std::max( fabs(tnew), fabs(t) ))) // check if close enough to be done
				return true;
			step = sgn(step) * d;
		}

		rk45::step( t, t+step, Yold, Ynew );
		// Check if that was the last step we needed to take
		if (integrationFlag == 0)
			return true;


		if ( epsilon > maxTolerance) {
			// Need to step back to current value of t and try again
			step *= 0.5;
			badStep++;
			if (badStep > maxIterations) // Give up if too many bad steps
				return false;
		} else {
			// step was good
			t = t+step;
			for (int i=0; i<neqns; i++)
				Yold[i] = Ynew[i];

			if ( epsilon < minTolerance )
				step *= 2;
		}
	}
	return true;
}

} /* namespace solvers */
