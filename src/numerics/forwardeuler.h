/*
 * forwardeuler.h
 *
 *  Created on: Nov 13, 2014
 *      Author: tony
 */

#ifndef FORWARDEULER_H_
#define FORWARDEULER_H_

#include "solver.h"

namespace solvers {


class forward_euler: public solver {
protected:
	double* dYdt;
public:
	forward_euler(cellmodel::Cells* c): solver(c){
		dYdt = new double[ neqns ];
	}
	virtual ~forward_euler(){
		delete[] dYdt;
	}
	virtual bool step( const double told, const double tnew, double* Yold, double* Ynew );
};

} /* namespace solvers */

#endif /* FORWARDEULER_H_ */
