/*
 * rk45.h
 *
 *  Created on: Nov 14, 2014
 *      Author: tony
 */

#ifndef RK45_H_
#define RK45_H_

#include "rk2.h"

namespace solvers {

class rk45: public rk2 {
protected:
	double epsilon;
	double *k2, *k3, *k4, *k5, *k6, *Y4;
public:
	rk45( cellmodel::Cells* c ):rk2( c ), epsilon(0){
		k2 = new double[ neqns ];
		k3 = new double[ neqns ];
		k4 = new double[ neqns ];
		k5 = new double[ neqns ];
		k6 = new double[ neqns ];
		Y4 = new double[ neqns ];
	}
	virtual ~rk45(){
		delete[] k2;
		delete[] k3;
		delete[] k4;
		delete[] k5;
		delete[] k6;
		delete[] Y4;
	}
	virtual bool step( const double told, const double tnew, double* Yold, double* Ynew );
private:
	static constexpr double c20 = 1.0/4.0, c21 = 1.0/4.0;
	static constexpr double c30 = 0.375, c31 = 0.09375, c32 = 0.28125;
	static constexpr double c40 = 12.0/13.0, c41 = 1932.0/2197.0, c42 = -7200.0/2197.0, c43 = 7296.0/2197.0;
	static constexpr double c51 = 439.0/216.0;
	static constexpr double c53 = 3680.0/513.0;
	static constexpr double c52 = -8.0, c54 = -845.0/4104.0;
	static constexpr double c60 = 0.5, c61 = -8.0/27.0, c62 = 2, c63 = -3544.0/2565.0, c64 = 1859.0/4104.0;
	static constexpr double c65 = -0.275;
	static constexpr double a1 = 25.0/216.0, a3 = 1408.0/2565.0, a4 = 2197.0/4104.0, a5 = -0.2;
	static constexpr double b1 = 16.0/135.0, b3 = 6656.0/12825.0, b4 = 28561.0/56430.0, b5= -0.18, b6 = 2.0/55.0;
};

} /* namespace solvers */

#endif /* RK45_H_ */
