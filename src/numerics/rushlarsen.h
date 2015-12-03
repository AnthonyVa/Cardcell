/*
 * rushlarsen.h
 *
 *  Created on: Dec 4, 2014
 *      Author: Anthony Varghese
 */

#ifndef RUSHLARSEN_H_
#define RUSHLARSEN_H_

#include "solver.h"

namespace solvers {

class rushlarsen: public solver {
protected:
	unsigned int MarkovIterations;
	double MarkovIterationError;
public:
	rushlarsen();
	rushlarsen( cellmodel::Cells* c ): solver( c ),
			MarkovIterations(0),
			MarkovIterationError(0){
	}
	virtual ~rushlarsen(){}
	virtual bool step( const double told, const double tnew, double* Yold, double* Ynew ){
		return cells->cellDERushLarsenUpdate(told, tnew);
	}
	// Rush Larsen solves for the new state in place - without using derivatives
	virtual bool isInplace(){
		return true;
	}

};

} /* namespace solvers */

#endif /* RUSHLARSEN_H_ */
