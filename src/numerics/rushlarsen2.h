/*
 * rushlarsen2.h
 *  Second-order Rush-Larsen technique as used in:
 *  "A Second-Order Algorithm for Solving Dynamic Cell Membrane Equations"
 *  Sundnes, J., Artebrant, R. ; Skavhaug, O. ; Tveito, A.
 *  Dept. of Inf., Univ. of Oslo, Oslo, Norway
 *  IEEE Trans Biomed Eng. 2009 Oct;56(10):2546-8. doi: 10.1109/TBME.2009.2014739. Epub 2009 Feb 20.
 *
 *  Created on: Dec 4, 2014
 *      Author: Anthony Varghese
 */

#ifndef RUSHLARSEN2_H_
#define RUSHLARSEN2_H_

#include <numerics/rushlarsen.h>

namespace solvers {

class rushlarsen2: public rushlarsen {
public:
	rushlarsen2();
	rushlarsen2( cellmodel::Cells* c ): rushlarsen( c ){
	}

	virtual ~rushlarsen2();

	virtual bool step( const double told, const double tnew, double* Yold, double* Ynew ){
		return cells->cellDERushLarsen2Update(told, tnew);
	}
	// Rush Larsen 2 solves for the new state in place - without using derivatives
	virtual bool isInplace(){
		return true;
	}
};

} /* namespace solvers */

#endif /* RUSHLARSEN2_H_ */
