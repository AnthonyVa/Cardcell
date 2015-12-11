/*
 * Cells.h
 *
 *  Created on: Dec 6, 2014
 *      Author: tony
 */

#ifndef CELLS_H_
#define CELLS_H_

#include <fstream>
using std::ofstream;
#include <cell/Moreno2011Cell.h>

class simulation;
enum Blockers{
	None,
	Lidocaine,
	Flecainide
};


namespace cellmodel {

class Cells {
public:
	Cells(const unsigned int n, const simulation* s, const double* par );
	void setSimulation( const simulation* sim );
	unsigned int getNcells() const {
		return ncells;
	}
	unsigned int getNstates() const {
		return nstates;
	}

	Moreno2011Cell getCell( const unsigned int n );
	void copyCell( const unsigned int i, const Moreno2011Cell& source);
	void copyCells( const Moreno2011Cell& source );
	void getStateFromCells ( double* Y );

	// Differential Equation Solver function
	// function that computes Cell Differential equations
	void celldiffeq(const double t, double* Y, double* dYdt);
	void SetInitialConditions();
	bool cellDERushLarsenUpdate(const double told, const double tnew);
	bool cellDERushLarsen2Update(const double told, const double tnew);
	void printV( ofstream* f , double t );
	void printNacell( ofstream* f, double t, unsigned int n );
	bool isPrintingNeeded();
	void printThreshold( unsigned int n, int cycle );
	void printCV( ofstream* f, int cycle );
	void resetPrintFlag();
	void ResetCells( double t );
	virtual ~Cells();

	// Things from simulation needed by Cell:
	int sim_type = 0;
	double waitTime = 0;
	double BCL = 0;
	double I_duration_S1=0, stimulus=0;
	Blockers blocker = None;
	bool HeartFailure = false;
	double dt = 0;
	double drug = 0;

	// Diffusion can be set using input file
	double Diffusion = 0.00154;		// CHANGE TO D / 2 IF HEART_FAILURE
	static constexpr double dx = 0.01;					//.01 //Can go from 0.01 - 0.02 cm
	static constexpr double Vthreshold = -75; // mV

private:
	Moreno2011Cell* cells;
	// Simulation parameters
	const simulation* sim;

	const unsigned int ncells;
	unsigned int nstates; // number of states per cell
	unsigned int nCellsPaced;

	void saveStateIntoCells( double* Y );
	void getRatesFromCells ( double* dYdt );
	void Calculate_Update( const unsigned int n );
	void Calculate_RushLarsenUpdate( const unsigned int n );
	void Calculate_RushLarsen2Update( const double told, const double tnew, const unsigned int n, const int step );
	void setStim(const double t);

	// Used to map big state vector into individual cells
	double** Yold2D = nullptr;
	double** Ynew2D = nullptr;

	// Used for 2nd order Rush Larsen
	double* Vhalf;

	static constexpr double time_diff_tol = 1e-6;
};

} /* namespace cellmodel */

#endif /* CELLS_H_ */
