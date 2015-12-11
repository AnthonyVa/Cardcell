/*
 * Cells.cpp
 *
 *  Created on: Dec 6, 2014
 *      Author: tony
 */
#include <omp.h>

#include <cell/Cells.h>
#include <sim/simulation.h>
#include <cmath>

#include <stdexcept>

#include <iostream>
using std::cerr;
using std::endl;

namespace cellmodel {

Cells::Cells( const unsigned int n, const simulation* s , const double* par) : ncells(n) {
	cells = new Moreno2011Cell[ncells];
	nstates = cells[0].nstates;
	setSimulation( s );
	nCellsPaced = s->nCellPaced;
	Diffusion = par[2];
	for (unsigned int i = 0; i < ncells; i++){
		cells[i].setXpos(i);
		cells[i].setParams( par );
		cells[i].setEnvironment( this );
	}
	Yold2D = new double*[nstates];
	Ynew2D = new double*[nstates];

	// Used for 2nd order Rush Larsen only
	//if ( sim->solve_type == rl2 )
	Vhalf = new double[ncells];
}

Cells::~Cells() {
	delete[] cells;
	delete[] Yold2D;
	delete[] Ynew2D;
	delete[] Vhalf;
}

void Cells::setSimulation( const simulation* s ){
	sim = s;

	sim_type = sim->sim_type;
	waitTime = sim->waitTime;
	BCL      = sim->BCL;
	I_duration_S1 = sim->I_duration_S1;
	stimulus = sim->stimulus;
	blocker  = sim->blocker;
	HeartFailure = sim->HeartFailure;
	dt		 = sim->dt;
	drug	 = sim->drug;
	// This is needed here to initialize cell parameters with env
	//  variables like drug conc, etc.
	// This code is duplicated in constructor - only one should
	//  abs. necessary but it is possible that this - setSim - may
	//  not get called.
	for (unsigned int i = 0; i < ncells; i++)
		cells[i].setEnvironment( this );
}


Moreno2011Cell Cells::getCell( const unsigned int n ){
	return cells[n];
}

void Cells::copyCell( const unsigned int i, const Moreno2011Cell& source){
	cells[i].Copy( source );
}

void Cells::copyCells( const Moreno2011Cell& source ){
	for (unsigned int i = 0; i < ncells; i++)
		copyCell( i , source);// References are messed up after reading structs from a file, so copy
}



void Cells::celldiffeq(const double time, double* Y, double* dYdt){
	for (unsigned int i=0; i<ncells; i++)
		cells[i].setI_stimulus( 0.0 );

	// cells to be paced
	if (sim_type == 1 || sim_type == 2) {
		const double cycle_time = fmod ( time, BCL );
		if (time > waitTime && cycle_time <= I_duration_S1)
			for (unsigned int i=0; i<ncells && i<nCellsPaced; i++)
				cells[i].setI_stimulus( stimulus ); // Apply stimulus to paced cells only
	}

	// Update the state in the cells data structure - don't really need to do this until we need to print things
	saveStateIntoCells( Y );

	try {
#pragma omp parallel for
		for (unsigned int i=0; i<ncells; i++) {
			// Compute membrane currents and intracellular eqns
			cells[i].Calculate_All(time);
			// Compute inter-cellular fluxes
			Calculate_Update( i );
		} // End of (parallel) for n
	} catch ( std::runtime_error& re ){
		cerr << re.what() << endl;
	} catch ( std::string& s) {
		cerr << s << endl;
	}
	getRatesFromCells( dYdt );
}

void Cells::setStim(const double t){
	for (unsigned int i=0; i<ncells; i++)
		cells[i].setI_stimulus( 0.0 );
	// cells to be paced
	if ( (sim_type == 1 || sim_type == 2) && (t > waitTime)) {
		const double cycle_time = fmod ( t, BCL );
		if ( (cycle_time > waitTime) &&
			 (cycle_time <= (waitTime+I_duration_S1) ))
			for (unsigned int i=0; i<ncells && i<nCellsPaced; i++)
				// Apply stimulus to paced cells only
				cells[i].setI_stimulus( stimulus );
	}
}
/*
 * First-order Rush Larsen
 */
bool Cells::cellDERushLarsenUpdate(const double told, const double tnew){
	setStim( tnew );
	try {
#pragma omp parallel for
		for (unsigned int i=0; i<ncells; i++) {
			// Compute membrane currents and intracellular eqns
			cells[i].computeRushLarsenStepState( told, tnew-told );
			// Compute inter-cellular fluxes
			Calculate_RushLarsenUpdate( i );
		} // End of (parallel) for n
		return true; // step was good
	} catch ( std::runtime_error& re ){
		cerr << re.what() << endl;
	} catch ( std::string& s) {
		cerr << s << endl;
	}
	return false;// step was not good
}

/*
 * Second-order Rush Larsen
 */
bool Cells::cellDERushLarsen2Update(const double told, const double tnew){
	const double dt = tnew-told;

	try {
		// Step to tn+1/2
		setStim( told + dt/2 );
		#pragma omp parallel for
		for (unsigned int i=0; i<ncells; i++) {
			// Compute membrane currents and intracellular eqns
			cells[i].computeRushLarsen2StepState( told, dt, 1, 0 );
			// Compute inter-cellular fluxes
			Calculate_RushLarsen2Update( told, tnew, i , 1);
		} // End of (parallel) for n

		// Step to tn+1
		setStim( tnew );
		#pragma omp parallel for
		for (unsigned int i=0; i<ncells; i++) {
			// Compute membrane currents and intracellular eqns
			cells[i].computeRushLarsen2StepState( told, dt, 2, Vhalf[i] );
			// Compute inter-cellular fluxes
			Calculate_RushLarsen2Update( told, tnew, i , 2);
		} // End of (parallel) for n

		return true; // step was good
	} catch ( std::runtime_error& re ){
		cerr << re.what() << endl;
	} catch ( std::string& s) {
		cerr << s << endl;
	}
	return false;// step was not good
}


void Cells::Calculate_Update( const unsigned int n ) {
	const double diffusion = Diffusion / (dx * dx);

	const double Vn   = cells[  n  ].getV();
	const double Vnm1 =    (n>0)     ? cells[n - 1].getV() : Vn;
	const double Vnp1 = (n<ncells-1) ? cells[n + 1].getV() : Vn;

	double Iaxnm1 = 0;
	double Iaxnp1 = 0;

	if ( ncells > 1 ){
		// At least 2 cells
		if ( n == 0 )
			Iaxnp1 = -(Vnp1 - Vn) * diffusion; // first cell
		else if ( n >= 1 && n < (ncells - 1)){
			Iaxnm1 = -(Vnm1 - Vn) * diffusion; // middle cells
			Iaxnp1 = -(Vnp1 - Vn) * diffusion; // middle cells
		} else if ( n == ncells - 1 )
			Iaxnm1 = -(Vnm1 - Vn) * diffusion; // last cell
	}
	cells[n].updatedVdt(Iaxnm1, Iaxnp1);
}


void Cells::Calculate_RushLarsenUpdate( const unsigned int n ) {
	const double diffusion = Diffusion / (dx * dx);

	const double Vn   = cells[  n  ].getV();
	const double Vnm1 =    (n>0)     ? cells[n - 1].getV() : Vn;
	const double Vnp1 = (n<ncells-1) ? cells[n + 1].getV() : Vn;

	double Iaxnm1 = 0;
	double Iaxnp1 = 0;

	if ( ncells > 1 ){
		// At least 2 cells
		if ( n == 0 )
			Iaxnp1 = -(Vnp1 - Vn) * diffusion; // first cell
		else if ( n >= 1 && n < (ncells - 1)){
			Iaxnm1 = -(Vnm1 - Vn) * diffusion; // middle cells
			Iaxnp1 = -(Vnp1 - Vn) * diffusion; // middle cells
		} else if ( n == ncells - 1 )
			Iaxnm1 = -(Vnm1 - Vn) * diffusion; // last cell
	}
	cells[n].updatedVdt(Iaxnm1, Iaxnp1);
	cells[n].fe_updateV( dt );
}

void Cells::Calculate_RushLarsen2Update( const double told, const double tnew,
										 const unsigned int n, const int step ) {
	const double diffusion = Diffusion / (dx * dx);
	const double delta = 1e-7; // for finite differencing
	const double Vn = cells[  n  ].getV();
	const double Vx = Vn + delta; // for finite-differencing
	const double dt = tnew-told;

	Calculate_Update( n ); // computes dV/dt for this cell
	const double a = cells[n].getdVdt();
	double b = 0;

	if (step == 1) {
		const double thalf = told + dt/2.0;
		// Finite difference to compute d(I_total)/dV
		const double Itotal	 = cells[n].getI_total(); // only true for 1 cell case
		const double Itotald = cells[n].CalculateCurrents( thalf, Vx);
		b = -(Itotald - Itotal)/delta;

		// n==1 --> Don't have to consider axial currents
		if ( ncells >= 2 ){
			// At least 2 cells
			if (( n == 0 ) || ( n == ncells-1 ))
				b = b - diffusion; // first or last cell
			else if ( n > 0 && n < (ncells - 1))
				b = b - 2 * diffusion; // middle cells
		}
		if (fabs(b) < delta)
			Vhalf[n] = Vn + a * dt/2;
		else
			Vhalf[n] = Vn + (a/b)*(exp(b*dt/2) - 1);
	}
	if (step == 2 ){
		// Finite difference to compute d(I_total)/dV
		const double Itotal	 = cells[n].getI_total(); // only true for 1 cell case
		const double Itotald = cells[n].CalculateCurrents( tnew, Vx);
		b = -(Itotald - Itotal)/delta;

		// n==1 --> Don't have to consider axial currents
		if ( ncells >= 2 ){
			// At least 2 cells
			if (( n == 0 ) || ( n == ncells-1 ))
				b = b - diffusion; // first or last cell
			else if ( n > 0 && n < (ncells - 1))
				b = b - 2 * diffusion; // middle cells
		}
		double Vnew = Vn;
		if (fabs(b) < delta)
			Vnew = Vn + a * dt;
		else
			Vnew = Vn + (a/b)*(exp(b*dt) - 1);
		cells[n].setV( Vnew );
	}
}

// loading in initial conditions
void Cells::SetInitialConditions(){
	for (unsigned int i = 0; i < ncells; i++)
		cells[i].SetInitialConditions();
}

void Cells::ResetCells( double t ){
	for (unsigned int i=0; i<ncells; i++)
		cells[i].Calculate_Reset(t);
}

void Cells::saveStateIntoCells( double* Y ){
	const int nstates = cells[0].nstates; // assumes all cells have same number of states

	// Update the state in the cells data structure
	/*
	 * First, set up a "2-d" array with nstates rows and ncells columns and set up
	 * pointers into the 1-d array Y
	 */
	for (int i=0; i<nstates; i++)
		Ynew2D[i] = &(Y[i*ncells]);

	/*
	 * Next, get a slice corresponding to the states of a single cell into the singlecell array
	 * and save each slice into a single cell
	 */
	double* singlecell = new double[nstates];
	for (unsigned int i=0; i<ncells; i++ ){
		for (int s=0; s<nstates; s++)
			singlecell[s] = Ynew2D[s][i];
		cells[i].setState( singlecell );
	}
	// Since nstates is around 55, these are not big arrays - creating and deleting should be fast
	delete[] singlecell;
}

void Cells::getStateFromCells ( double* Y ){
	const int nstates = cells[0].nstates; // assumes all cells have same number of states

	// Copy state information from the cells to an array
	// First, set up a "2-d" array with nstates rows and ncells columns and set up
	//  pointers into the 1-d array Y
	for (int i=0; i<nstates; i++)
		Yold2D[i] = &(Y[i*ncells]);

	// Next, get a slice corresponding to the states of a single cell into the singlecell array
	//  and save each slice into a single cell
	double* singlecell = new double[nstates];
	for (unsigned int i=0; i<ncells; i++ ){
		cells[i].getState( singlecell );
		for (int s=0; s<nstates; s++)
			Yold2D[s][i] = singlecell[s];
	}

	// Since nstates is around 51, these are not big arrays - creating and deleting should be fast
	delete[] singlecell;
}


void Cells::getRatesFromCells ( /*const*/ double* dYdt ){
	const int nstates = cells[0].nstates; // assumes all cells have same number of states

	// Copy state information from the cells to an array
	/*
	 * First, set up a "2-d" array with nstates rows and ncells columns and set up
	 * pointers into the 1-d array Y
	 */
	double** dYdt2D = new double*[nstates];
	for (int i=0; i<nstates; i++)
		dYdt2D[i] = &(dYdt[i*ncells]);

	/*
	 * Next, get a slice corresponding to the states of a single cell into the singlecell array
	 * and save each slice into a single cell
	 */
	double* singlecell = new double[nstates];
	for (unsigned int i=0; i<ncells; i++ ){
		cells[i].getRates( singlecell );
		for (int s=0; s<nstates; s++)
			dYdt2D[s][i] = singlecell[s];
	}
	// Since nstates is around 55, these are not big arrays - creating and deleting should be fast
	delete[] singlecell;
	delete[] dYdt2D;
}

void Cells::printV( ofstream* f , double t ){
	*f << t << ", ";
	for (unsigned int i=0; i<ncells; i++)
		*f << cells[i].getV() << ", ";
	*f << std::endl;
}

void Cells::printNacell( ofstream* f, double t, unsigned int n ){
	*f << t << ", "
			<< cells[n].getNaClosed()	<< ", "
			<< cells[n].getNaIC23()		<< ", " << cells[n].getNaIF() << ", "
			<< cells[n].getNaOS()		<< ", " << cells[n].getNaO()  << ", "
			<< cells[n].getNaDClosed()	<< ", "
			<< cells[n].getNaDIC23()	<< ", "	<< cells[n].getNaDIF()	<< ", "
			<< cells[n].getNaDOS()		<< ", "	<< cells[n].getNaDO()	<< ", "
			<< cells[n].getNaDIM1()		<< ", "
			<< cells[n].getNaD_Closed()	<< ", "
			<< cells[n].getNaD_IC23()	<< ", "	<< cells[n].getNaD_IF()	<< ", "
			<< cells[n].getNaD_OS()		<< ", "	<< cells[n].getNaD_O()	<< ", "
			<< cells[n].getNaD_IM1()	<< std::endl;
}

// check if all cell flag2's are true
bool Cells::isPrintingNeeded(){
	for (unsigned int i=0; i<ncells; i++ ) {
		if ( !cells[i].flag2 )
			return false; // no need to check any more cells, get out of for loop
	}
	return true;
}



void Cells::printThreshold( unsigned int n, int cycle ){
	if (cells[50].V_thr < Vthreshold) {
		std::cout << "Threshold is less than " << Vthreshold
				<< " at Cell 50 beat number: " << cycle
				<< ", " << "Drug Concentration is:  "
				<< sim->drug << std::endl;
	}

}


void Cells::printCV( ofstream* f, int cycle ){

	const double conversion  = 1000; // cm/sec --> 1000 msec cm / sec
	const double cell_length = 0.01; // 100 um = 0.1 mm = 0.01 cm

	// Compute conduction velocity only for chains of length > 24
	if ( ncells >25 ){
		// 3 ways to measure conduction velocity:
		// assumes wave is traveling from left to right - low indices to high

		// 1. adjacent cells in middle
		const unsigned int middle = ncells/2;
		const double time_diff1 = cells[ middle ].t_thr - cells[ middle-1 ].t_thr;
		if ( fabs(time_diff1) < time_diff_tol || fabs(time_diff1) > BCL )
			return;
		const double length1 = cell_length;
		const double cv1 = conversion * length1/time_diff1; // cm/msec

		// 2. middle to 3/4 cell
		const unsigned int three4th = 3*ncells/4;
		const double time_diff2 = cells[ three4th ].t_thr - cells[ middle ].t_thr;
		if ( fabs(time_diff2) < time_diff_tol )
			return;
		const double length2 = (three4th - middle) * cell_length;
		const double cv2 = conversion * length2/time_diff2; // cm/msec

		// 3. 3/4 cell
		const double time_diff3 = cells[ three4th ].t_thr - cells[ three4th-1 ].t_thr;
		if ( fabs(time_diff3) < time_diff_tol )
			return;
		const double length3 = cell_length;
		const double cv3 = conversion * length3/time_diff3; // cm/msec

		*f << cycle << ", " << cv1 << ", " << cv2 << ", " << cv3 << std::endl;
	}
}


void Cells::resetPrintFlag(){
	for (unsigned int i=0; i<ncells; i++ )
				cells[i].flag2 = false;
}

} /* namespace cellmodel */
