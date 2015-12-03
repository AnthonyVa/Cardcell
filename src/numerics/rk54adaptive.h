/*
 * rk54adaptive.h
 *
 *  Created on: Nov 14, 2014
 *      Author: W1084063
 */

#ifndef RK54ADAPTIVE_H_
#define RK54ADAPTIVE_H_

#include "rk45.h"

namespace solvers {

class rk54adaptive: public rk45 {
protected:
	double minTolerance, maxTolerance, minStep, maxStep;
public:
	rk54adaptive( cellmodel::Cells* c,
				const double epsmin,	const double epsmax,
				const double hmin,		const double hmax) :	rk45( c ),
				minTolerance( epsmin ),	maxTolerance( epsmax ),
				minStep( hmin ),		maxStep( hmax ){
	}
	virtual ~rk54adaptive(){
	}
	virtual bool step( const double told, const double tnew, double* Yold, double* Ynew );
private:
	static const int maxIterations = 1000;
	static constexpr double delta = 1e-6; // units are milliseconds
};

} /* namespace solvers */

#endif /* RK54ADAPTIVE_H_ */
