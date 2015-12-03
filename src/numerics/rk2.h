/*
 * rk2.h
 *
 *  Created on: Nov 13, 2014
 *      Author: tony
 */

#ifndef RK2_H_
#define RK2_H_

#include "forwardeuler.h"

namespace solvers {

class rk2: public forward_euler {
protected:
	double* k1;
	double* Yhalf;
public:
	rk2(cellmodel::Cells* c):forward_euler( c ){
		k1    = new double[ neqns ];
		Yhalf = new double[ neqns ];
	}
	virtual ~rk2(){
		delete[] k1;
		delete[] Yhalf;
	}
	virtual bool step( const double told, const double tnew, double* Yold, double* Ynew );
};

} /* namespace solvers */

#endif /* RK2_H_ */
