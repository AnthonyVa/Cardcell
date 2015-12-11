/*
 * simulation.cpp
 *
 *  Created on: Nov 23, 2014
 *      Author: tony
 */

#include <iostream>

using std::cin;
using std::cout;
using std::endl;
using std::flush;
using std::string;

#include <sim/simulation.h>




simulation::simulation(const unsigned int type, const unsigned int n,
		const double* par, const string folder, const string fpre,
		const unsigned int outlev) :
				sim_type( type ),
				outlevel(outlev),
				nCell( n ),
				folderName(folder){
	mycells = new cellmodel::Cells( nCell, this , par );
	if (nCell==1) {
		nCellOutputs = 1;
		nCellsToOutput = new unsigned int[nCellOutputs];
		nCellsToOutput[0] = 0;
		files   = new filesStruct(folder, fpre, 1 );
	} else {
		nCellOutputs = 5;
		nCellsToOutput = new unsigned int[nCellOutputs];
		nCellsToOutput[0] = 0;
		nCellsToOutput[1] = nCell/(nCellOutputs-1);
		for (unsigned int i=2; i<nCellOutputs-1; i++)
			nCellsToOutput[i] = nCellsToOutput[i-1] + nCellsToOutput[1];
		nCellsToOutput[nCellOutputs-1] = nCell-1;
		files   = new filesStruct(folder, fpre, nCellOutputs /* # of cell files */);
	}
	// The next line assumes all cells have the same number of states
	const unsigned int nstates = mycells->getNstates();
	neqns = nstates * nCell;
	Yold	= new double[neqns];
	Ynew	= new double[neqns];
}

simulation::~simulation(){
	// Deallocate work space
	delete[] Yold;
	delete[] Ynew;

	delete[] nCellsToOutput;
	delete files;
	delete mycells;
}

void simulation::SetInitialConditions( ){
	// Load in initial conditions for each cell in the array
	mycells->SetInitialConditions();

	switch (sim_type){
	case 0:
		// Sim type 0 is for initial conditions set in code - no file read.
		break;
	case 1:
		// Read in file containing state of 1 cell
		// Sim type 1 typically uses the state of a single cell to set
		//  the states of all the cells in a cable
		// We need to set the initial conditions first - done above - and then
		//  reset the values using the data from the input file
		// This is because we use references in the Cell class/struct
		// Reading references from the input file does not work.
		readsinglecell( );
		break;
	case 2:
		// Sim type 2 starts up from  last cable state - read in the
		// states of all cells from a file
		readcells( );
		break;
	default:
		break;
	}

	// Get a complete copy of all the state variables of all the cells
	mycells->getStateFromCells( Yold );
}



void simulation::printcellAux(ofstream& of, const double t, Moreno2011Cell& c) {
	of << t << ", "
			<< c.getV() << ", " << c.getI_K1() << ", " << c.getI_Na() << ", "
			<< c.getI_axnm1() << ", " << c.getI_axnp1() << ", " << c.getI_axial()
			<< endl;
}

void simulation::printcellEAD(ofstream& of, Moreno2011Cell& c, Moreno2011Cell& c2, int r, int p, double S1) {
	double th_diff = 0;
	if (r == 0)
		th_diff = 1000 * 0.01 / (c.t_thr - c2.t_thr);
	else
		th_diff = 1000 * 0.01 / (c2.t_thr - c.t_thr); // special case for Cell[0]

	of << c.t_min << ", " << c.V_min << ", " << c.t_thr << ", " << c.V_thr
			<< ", " << c.t_max << ", " << c.V_max << ", " << c.t_90 << ", "
			<< c.V_90 << ", " << c.t_EAD << ", " << c.V_EAD << ", " << c.t_EAD2
			<< ", " << c.V_EAD2 << ", " << c.L_EAD << ", " << c.APD_90 << ", "
			<< c.DI << ", " << p << ", " << S1 << " , " << th_diff << ", "
			<< c.peak_slope << endl;
}

/*
 * savecells
 */
void simulation::savecells() {
	// if n == 1 --> save a single cell state
	//    n > 1  --> save multiple cells state
	const unsigned int ncells = mycells->getNcells();
	string f = (ncells==1)	? (folderName + "/" + cellstateFileName)
									: (folderName + "/" + stateFileName);
	if (outlevel > 1)
		cout << "    Writing to: " << f << ". " << flush;

	FILE *fp = fopen(f.c_str(), "wb");
	int nw=0;
	for (unsigned int i = 0; i < ncells; i++) {
		Moreno2011Cell temp = mycells->getCell( i );
		nw += fwrite(temp.state, sizeof(Moreno2011Cell::state), 1, fp);
	}
	fclose(fp);

	if (outlevel > 1)
		cout << " Wrote state of " << nw << " cell(s)" << std::endl << flush;
}

/*
 * readcells
 */
void simulation::readcells() {
	string f = folderName + "/" + cellstateFileName;
	const unsigned int ncells = mycells->getNcells();
	if (outlevel > 1)
		cout << "    Reading data for " << ncells << " cells from " << f << ".";

	FILE *fp = fopen(f.c_str(), "r");
	int nr = 0;
	unsigned int ntotal = 0;
	// References are messed up after reading structs from a file
	for (unsigned int i = 0; i < ncells; i++) {
		Moreno2011Cell temp;
		nr = fread(temp.state, sizeof(Moreno2011Cell::state), 1, fp);
		if (nr <= 0)
			break;
		mycells->copyCell( i, temp);
		ntotal++;
	}
	fclose(fp);

	if (outlevel > 1){
		if ( ntotal < ncells )
			cout << " --> Read " << nr << " cells. ";
		cout << " Done reading." << endl << flush;
	}
}

/*
 * readsinglecell
 */
void simulation::readsinglecell() {
	string f = folderName + "/" + cellstateFileName;
	const unsigned int ncells = mycells->getNcells();
	if (outlevel > 1)
		cout << "    Reading data for " << ncells << " cells from " << f << ".";

	FILE *fp = fopen(f.c_str(), "r");
	Moreno2011Cell temp;
	int nr = fread(temp.state, sizeof(Moreno2011Cell::state), 1, fp);
	switch (nr){
	case 0:
		cout << " --> Could not read cell state!!! ";
		break;
	case 1:
		if (outlevel > 1)
			cout << " --> Read single cell state. ";
		mycells->copyCells( temp );
		break;
	default:
		cout << " --> fread returned: " << nr << ", cell state not changed.";
		break;
	}
	fclose(fp);

	if (outlevel > 1)
		cout << " Done setting cell states from file." << endl << flush;
}



void simulation::printEverything(const double time, const int cycle,
					 const double interval, const int S2_cycle,
					 int& counter ){
	counter++;
	if ((counter % counterupdates) == 0) {
		if (outlevel > 1)
			cout << ".";
		// Print auxilliary variables
		if ((counter % (printtimeupdate*counterupdates)) == 0){
			//This file will have voltages against time for each cell
			mycells->printV( files->result, time );
			if (printMarkovIts)
				printMarkovIterations( time );

			//This file prints The Markov States for cell 50
			if (cycle > printMarkovStatesAfterCycle) {
				mycells->printNacell( files->Na, time, nCell/2 );
			}

			if (outlevel > 1)
				cout << "t=" << time << "msec" << endl;
			for (unsigned int i=0; i<nCellOutputs; i++){
				Moreno2011Cell  c = mycells->getCell( nCellsToOutput[i] );
				printcellAux( *(files->cellAuxFs[ i ]), time, c );
			}

			/* Used to be:
			Cell  c0 = mycells->getCell( n0p00 );
			printcellAux( *(files->First), time, c0 );
			if (n0p25 > n0p00){
				Cell c25 = mycells->getCell( n0p25 );
				printcellAux( *(files->Quarter), time, c25 );
				if (n0p50 > n0p25){
					Cell c50 = mycells->getCell( n0p50 );
					printcellAux( *(files->Midpoint), time, c50 );
					if (n0p75 > n0p50){
						Cell c75 = mycells->getCell( n0p75 );
						printcellAux( *(files->Three4ths), time, c75 );
						if (n1p00 > n0p75){
							Cell c99 = mycells->getCell( n1p00 );
							printcellAux( *(files->Last), time, c99 );
						}
					}
				}
			}
			*/
		}

		if ( (counter % (printCycleUpdate*printtimeupdate*counterupdates)) == 0){
			if (outlevel > 1)
				switch (sim_type){
				case 0:
					cout << "Non-paced simulation, t=" << time  << " end time= ";
					break;
				case 1:
				case 2:
					cout << "Paced simulation, cycle: " << cycle << " of " << MaxNumberOfBeats
							<< ", t=" << time  << ", this interval= ";
					break;
				default:
					break;
				}
			counter = 0;
			if (outlevel > 1)
				cout  << interval << ", [drug] is " << drug << endl;
		}
		cout << flush;
	}


	if ( nCell > 1 ){
		//This file is for individual cells that calculates all of the min, max, V_90 parameters.
		bool isPrintingNeeded = mycells->isPrintingNeeded();
		if (isPrintingNeeded) {
			mycells->printThreshold( 50, cycle );

			for (unsigned int i=0; i<nCellOutputs; i++){
				unsigned int cellout  = nCellsToOutput[i];
				unsigned int cellnext = 1;
				int flag = 0;
				if (nCellsToOutput[i] == 0 ){
					// For the 0th cell, the adjacent one is 1
					cellnext = nCellsToOutput[i] + 1;
					flag = 1;
				} else {
					// For the other cells, the adjacent one is the one before it
					cellnext =  nCellsToOutput[i] - 1;
				}

				Moreno2011Cell  c    = mycells->getCell( cellout );
				Moreno2011Cell  cadj = mycells->getCell( cellnext );

				printcellEAD( *(files->EADFs[ i ]),  c,  cadj, flag, cellout,  cycle );
			}

			/* Used to be
			Cell  c0 = mycells->getCell( n0p00 );
			Cell  c1 = mycells->getCell( n0p00+1 );
			printcell( *(files->EAD00),  c0,  c1, 1,     0,   cycle );
			Cell c25 = mycells->getCell( n0p25 );
			Cell c24 = mycells->getCell( n0p25-1 );
			printcell( *(files->EAD25), c25, c24, 0,    25,   cycle );
			Cell c50 = mycells->getCell( n0p50 );
			Cell c49 = mycells->getCell( n0p50-1 );
			printcell( *(files->EAD50), c50, c49, 0,    50,   cycle );
			Cell c75 = mycells->getCell( n0p75 );
			Cell c74 = mycells->getCell( n0p75-1 );
			printcell( *(files->EAD75), c75, c74, 0,    75,   cycle );
			Cell c99 = mycells->getCell( n1p00 );
			Cell c98 = mycells->getCell( n1p00-1 );
			printcell( *(files->EAD99), c99, c98, 0, nCell-1, cycle );
			*/

			mycells->resetPrintFlag();
			isPrintingNeeded = false;
		} // End of individual result files
	}
}

bool simulation::step( const double told, const double tnew ){
	const bool success = solve->step( told, tnew, Yold, Ynew );
	if (success) {
		if (solve->isInplace()) {
			// Get a copy of all the state variables of all the cells
			mycells->getStateFromCells( Yold );
		} else {
			// Update the old state to have the new values
			for (unsigned int i=0; i<neqns; i++)
				Yold[i] = Ynew[i];
		}
	}
	return success;
}

void simulation::printCV(const double time, const int cycle){
	if ( nCell > 1 ){
		mycells->printCV( files->CV, cycle );
	}
}
void simulation::printMarkovIterations( const double t ){
	const unsigned int ncells = mycells->getNcells();
	// find middle cell
	const unsigned int middlecell = ncells/2;
	Moreno2011Cell  c = mycells->getCell( middlecell);
//MarkovIterations should be in the rushlarsen object
	*(files->MarkovIterations) << t <<", " << c.MarkovIterations
						       << ", " << c.MarkovIterationError << endl;
	// reset markov iteration count
}
