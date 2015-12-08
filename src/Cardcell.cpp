//============================================================================
// Name        : Cardcell.cpp
// Author      : Anthony Varghese
// Version     :
// Copyright   : All rights reserved by author(s)
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
using std::cout;
using std::endl;
using std::flush;

#include <sim/simulation.h>

/********************************************************************/
// File I/O Functions
const char* inputfile = "inputs/singlecellS.inputparams";
void readInputParameters( const char * inputfile, simulation**& sim,
						  unsigned int& nsims );
double getRealTime();
/********************************************************************/


int main() {
	cout << "Simple cardiac cell!" << endl; // prints Simple cardiac cell!

	//Initializations
	unsigned int NumberOfSimulations = 0;
	simulation** sims = nullptr; // array of pointers to simulation
	readInputParameters( inputfile, sims, NumberOfSimulations );
	if (sims == nullptr) {
		cout << "Quitting program - Simulations not set up!" << endl;
		return -1;
	}

	for (unsigned int nsim = 0; nsim < NumberOfSimulations; nsim++) {
		simulation* s = sims[ nsim ];

		const double start_time = getRealTime();

		// Load in initial conditions for each cell in the array
		s->SetInitialConditions();

		const unsigned int S2_cycle = 0; // not really used currently
		const unsigned int LastCycle = s->MaxNumberOfBeats;
		double t = 0;
		int stepCounter = -1; // Get's incremented in printEverything
		bool stop = false;
		s->printEverything(0, 0, s->waitTime, 0, stepCounter );

		//	 	 Begin Time Loop Here
		for (unsigned int S1_cycle = 1; S1_cycle <= LastCycle && !stop;
							S1_cycle++) {
			const double interval = ( (s->sim_type==0) && (S1_cycle == 1) )
					? s->waitTime : s->BCL;
			const double dt = s->dt;
			const int sim_type = s->sim_type;

			for (double tTemp = 0; tTemp <= interval && !stop; tTemp += dt) {
				const double tnew = t + dt;
				const bool successfulstep =	s->step( t, tnew );
				if (successfulstep)
					t = tnew;

				// Check if we need to save the state of the cells
				if (tTemp >= (interval - dt))
					// Reached the end of the simulation
					if ( (sim_type==0) ||
						 ((sim_type==1 || sim_type==2) && (S1_cycle==LastCycle)) ){
						s->savecells();	// Single or multiple cells
						stop = true;
						continue;
					}
				s->printEverything(t, S1_cycle, interval, S2_cycle, stepCounter );
			} //End of For tTemp ==
			s->printCV( t, S1_cycle);
			//Reset all cell characteristics such as APD, DI, V90 etc. every beat
			s->mycells->ResetCells( t );
		} //End of for S1

		delete s; // to save memory
		const double end_time = getRealTime( );
		cout << "Timing Information: " << (end_time - start_time) << " seconds" << endl;
	} // sim loop

	delete[] sims;

	return 0;
}
