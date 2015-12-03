/*
 * solver.h
 *
 *  Created on: Nov 13, 2014
 *      Author: tony
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include <cell/Cells.h>

namespace solvers {

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

class solver {
protected:
	cellmodel::Cells* cells;
	int nstates; // assumes all cells have same number of states
	int neqns;   // total number of equations
public:
	solver(cellmodel::Cells* c):	cells( c ),
						nstates( c->getNstates() ),
						neqns( c->getNcells() * nstates ) {
	}
	virtual ~solver();
	// Some solvers (like Rush Larsen) solve for the new state in place - without using derivatives
	// Most solvers are not in place solvers
	virtual bool isInplace(){
		return false;
	}
	virtual bool step( const double told, const double tnew, double* Yold, double* Ynew ){ return false;}
};

} /* namespace solvers */

#endif /* SOLVER_H_ */
