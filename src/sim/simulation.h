/****************************************************************************
 *
 *  Based on code from Moreno et al. 2011:
 *
 * 				Global Variables for Human Cable Simulations
 * 				By: Jonathan D. Moreno
 * 				Colleen Clancy Laboratory
 * 				August 2007
 *
 *  Created on: Sep 30, 2014
 *      Author: Anthony Varghese
 *
 *****************************************************************************/


#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <sstream>
#include <fstream>
using std::string;
using std::stringstream;
using std::ofstream;


#include <cell/Cells.h>
// Numerical integrators
#include <numerics/rk45.h>
#include <numerics/rk54adaptive.h>
using namespace solvers;


struct filesStruct{
	ofstream* result;
	ofstream* Na;
	unsigned int ncellfiles;
	ofstream** cellAuxFs;

	ofstream* CV;
	// These files hold data for individual cells to compare V_min, V_max, V_90 etc. to compare between cells
	ofstream** EADFs;
	ofstream* MarkovIterations;

	filesStruct(const string folder, const string fpre, const unsigned int ncellf){
		string prefix = folder + "/" + fpre;
		result	 = new ofstream(( prefix + "V.txt").c_str() );
		Na		 = new ofstream(( prefix + "Na_middleCell.txt").c_str());

		ncellfiles = ncellf;
		cellAuxFs = new ofstream*[ncellfiles];
		for (unsigned int i=0; i<ncellfiles; i++){
			stringstream sstm;
			sstm << prefix << "cell" << i << ".txt";
			cellAuxFs[i] = new ofstream( (sstm.str()).c_str() );
		}

		CV	  = new ofstream(( prefix + "CV.txt").c_str());

		EADFs = new ofstream*[ncellfiles];
		for (unsigned int i=0; i<ncellfiles; i++){
			stringstream sstm;
			sstm << prefix << "cellEAD" << i << ".txt";
			EADFs[i] = new ofstream( (sstm.str()).c_str() );
		}

		MarkovIterations = new ofstream(( prefix + "MarkovIter.txt").c_str());
	}
	~filesStruct(){
		result->close();	delete result;
		Na->close();		delete Na;
		for (unsigned int i=0; i<ncellfiles; i++) {
			cellAuxFs[i]->close();
			delete cellAuxFs[i];
		}
		delete[] cellAuxFs;
		CV->close();		delete CV;

		for (unsigned int i=0; i<ncellfiles; i++) {
			EADFs[i]->close();
			delete EADFs[i];
		}
		delete[] EADFs;
		MarkovIterations->close(); delete MarkovIterations;
	}
};



class simulation{
public:
	// Constructor - the argument n is the number of cells
	simulation(const unsigned int type, const unsigned int n,
			const double* par, const string folder, const string fpre,
			const unsigned int outlev);
	void set_params( const double dt_in, const double BCL_in,
					 const unsigned int MaxCycles, const double WaitTime,
					 const double stim, const double s1dur,
					 const unsigned int ncellPace, const string blckr,
					 const double drg, const unsigned int counterUpd,
					 const unsigned int ptupd, const unsigned int pCycupd,
					 const unsigned int pMStAfter, const bool logMark){
		dt = dt_in;
		BCL = BCL_in;
		MaxNumberOfBeats = MaxCycles;
		waitTime = WaitTime;
		stimulus = stim;
		I_duration_S1 = s1dur;
		nCellPaced = ncellPace;
		if (blckr.compare("") == 0 )
		   blocker = None;
		else if (blckr.compare("Lidocaine") == 0 )
		   blocker = Lidocaine;
		else if (blckr.compare("Flecainide") == 0 )
		   blocker = Flecainide;
		drug = drg;
		counterupdates = counterUpd;
		printtimeupdate = ptupd;
		printCycleUpdate = pCycupd;
		printMarkovStatesAfterCycle = pMStAfter;
		printMarkovIts = logMark;
		mycells->setSimulation( this );
	}
	void setFolder( const string str ){
		folderName.append( str );
	}
	void setSolver( solver* const s){ solve = s; }
	void SetInitialConditions();
	virtual bool step( const double told, const double tnew );

	void printEverything(const double time, const int cycle,
						 const double interval, const int S2_cycle,
						 int& counter );
	void printCV(const double time, const int cycle);
	void printcellAux( ofstream& of, const double t, Moreno2011Cell& c );
	void printcellEAD( ofstream& of, Moreno2011Cell& c , Moreno2011Cell& c2, int r, int p, double S1);
	void savecells( );
	void readcells( );
	void readsinglecell( );
	void printMarkovIterations( const double t );
	virtual ~simulation();

	// Needed by numerical integrators to access mycells->diffeq
	cellmodel::Cells* mycells = nullptr;
	// Stimulus parameters
	unsigned int nCellPaced = 4;
	double I_duration_S1 = 0.2;
	double I_duration_S2 = 2.0;
	// threshold for a single cell - Lidocaine model:
	//                      BCL
	// duration   1000		500		350
	//	0.1ms      274		250		260
	//	0.2ms	   137		140		141
	//	0.5ms	    56		 56		 58
	//	1 ms	    29		 28		 30
	//	2 ms	    15		 15		 15
	//	5 ms	     6.5      6.6	  7.0
	//
	// Threshold for stimulating a single cell in chains of various sizes BCL == 1000 ms
	//			Number of cells in chain
	//				1		2		3		5		8		16
	// duration
	//	0.1ms      255	  497	  727	   987	 1072	  1090
	//	0.2ms	   140	  273	  405	   577    631	   638
	//	0.5ms	    55	  111	  165	   257	  298	   302
	//	1 ms	    29	   57	   85	   137    176	   181
	//	2 ms	    15	   30	   44	    73	  105	   114
	//	5 ms	     6.5   13.5    20		33	   51	    67

	// Threshold for stimulating cells in a chain of 100 cells BCL == 1000 ms
	//			Number of cells stimulated
	//				1		2		4		8
	// duration
	//	0.1ms     1100	  645	  435	  326
	//	0.2ms	   640	  370	  245	  181
	//	0.5ms	   302	  168	  105	   75
	//	1 ms	   185	   98	   58	   40
	//	2 ms	   115	   60	   34	   22
	//	5 ms	    67	   35	   19	   12
	double stimulus = 1.15 * -245.0; // 1.15 * -150; // normalized to Cm

	Blockers blocker = Lidocaine; // Flecainide;
	double drug= 0.0; // 20.0e-6;

	//const double waitTime = 600000;	// If sim_type ==0, use 10 minutes ==> 10 min * 60 sec/min ==> 600,000ms for wait time
	double waitTime = 1;	// If sim_type ==0, use 10 minutes ==> 10 min * 60 sec/min ==> 600,000ms for wait time
	//double waitTime = 0.1;		// if sim_type ==1, use 0.1
	double MaxNumberOfBeats = 100; // Number of beats or stimuli to simulate
	double BCL      = 1000.0;		// in milliseconds - time between stimuli


	// Numerical integration computation, progress updates, and file i/o
	double dt = 1.0e-6;  // Drug concentration in M
	int counterupdates   = 5000;					// print a . every so often
	int printtimeupdate  =  100 * counterupdates;	// print time every so often
	int printCycleUpdate =  10 * printtimeupdate;

	bool printMarkovIts = false;
	int printMarkovStatesAfterCycle = MaxNumberOfBeats - 2; // Print to Na.txt after this cycle number

	// Types of simulations:
	//  0 is absolute initials = ic set in code
	//  1 is for initial conditions read from single cell file
	//  2 is for initial conditions read from cable file
	unsigned int sim_type = 0;
	bool HeartFailure = false;	// HeartFailure == true --> parameters for Heart failure, else, regular parameters

private:
	unsigned int outlevel = 1;
	// Number of cells in array
	unsigned int nCell = 1;
	unsigned int nCellOutputs = 1;
	unsigned int* nCellsToOutput = nullptr;
//	unsigned int TimedUpdate = 1000; // Print output every 1000 msec - i.e. 1 second
	// File Names
	const string cellstateFileName  = "cellInitial_Conditions.dat"; // This file is used to store the initial hold of the 100 cell cable
	const string stateFileName		= "Initial_Conditions1.dat"; // This file is used to store the initial hold of the 100 cell cable
	const string stateFileName2		= "Initial_Conditions2.dat"; // This file is used to store the conditions after a certain number of paced beats

	string folderName;
	filesStruct* files;

	solver* solve = nullptr;
	unsigned int neqns;
	double* Yold = nullptr;
	double* Ynew = nullptr;
};

#endif
