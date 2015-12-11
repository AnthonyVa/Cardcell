/*
 * readparams.cpp
 *
 *  Created on: Dec 9, 2014
 *      Author: Anthony Varghese
 */

#include <iostream>
#include <fstream>
#include <cstdlib> // for system function
#include <json/json.h>
#include <numerics/solver.h>
#include <numerics/rk45.h>
#include <numerics/rk54adaptive.h>
#include <numerics/rushlarsen2.h>
using namespace solvers;

#include <cell/Cells.h>
#include <sim/simulation.h>

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::string;

/*
 * readInputParameters - read input parameters from json file into sim list
 */
void readInputParameters( const char * inputfile, simulation**& sims, unsigned int& nsims ){
	// Create an input file stream object to read in the
	//  text of the inputfile and then parse as a json file
	std::ifstream ifs(inputfile);
	string document((std::istreambuf_iterator<char>(ifs)),
						  std::istreambuf_iterator<char>());
	Json::Reader reader;
	Json::Value  root;
	bool parsedSuccessful = reader.parse( document, root );
	if (!parsedSuccessful){
		cerr << "Not able to parse input file!" << endl;
		return;
	}

	// At this point, the input has been parsed, we can start
	//  extracting all the data.
	nsims = root.get("NumberOfSimulations" , "0" ).asInt();
	if (nsims <= 0)
		return;

	std::string folder = root.get("OutputFolder" , "" ).asString();
	if (folder.compare( "" ) == 0)
		folder = "."; // Output will go into the current folder - project dir
	else {
		// Alternative to stuff below:
		// mkdir( str) - POSIX
		// or _mkdir( str ) in Win32 -- #include direct.h, stdlib.h
#if defined(_WIN32)
		string command = "mkdir " + folder;
#else
		string command = "mkdir -p " + folder;
#endif
		system(  command.c_str() );
	}
	cout << "Read input file; "  << nsims
			<< " simulations, output will go into: "<< folder << endl;

	sims  = new simulation*[nsims];

	const Json::Value simulations = root["Simulations"];
	const Json::Value cellparams  = root["Cells"];
	// Set cell parameters first - need a better way to set parameters.
	double params[3];
	params[0] = cellparams.get("GNa", "1.0").asDouble();
	params[1] = cellparams.get("GK1", "1.0").asDouble();
	params[2] = cellparams.get("Diffusion", "1.0").asDouble();
	for ( unsigned int index = 0; index < nsims; ++index ){
		// Iterates over the sequence elements.

	   int sim_type = simulations[index].get("type", "0" ).asInt();
	   int ncells = simulations[index].get("ncells", "1" ).asInt();
	   string fprefix = simulations[index].get("outputfileprefix", "").asString();

	   sims[index] = new simulation(sim_type, ncells, params, folder, fprefix);
	   simulation* s = sims[index];


		// Set up solver
	   string solvername = simulations[index].get("solver", "fe").asString();
	   solver* slvr = nullptr;
	   if (solvername.compare("fe") == 0){
			slvr = new solvers::forward_euler( s->mycells );
	   } else if (solvername.compare("rk2") == 0){
			slvr = new solvers::rk2( s->mycells );
	   } else if (solvername.compare("rk45") == 0){
			slvr = new solvers::rk45( s->mycells );
	   } else if (solvername.compare("rk45ad") == 0) {
		   slvr  = new solvers::rk54adaptive( s->mycells,
											1.0e-8,	// epsmin
											1.0e-4, // epsmax
											1.0e-6, // hmin
											1.0e-2  // hmax
		   );
	   } else if (solvername.compare("rl") == 0) {
		   slvr = new solvers::rushlarsen( s->mycells );
	   } else if (solvername.compare("rl2") == 0) {
		   slvr = new solvers::rushlarsen2( s->mycells );
	   }
	   s->setSolver( slvr );

	   double dt	 =	simulations[index].get("dt", 1e-6 ).asDouble();
	   double BCL	 =	simulations[index].get("BCL", 1000.0 ).asDouble();
	   int MaxCycles =	simulations[index].get("MaxCycles", 1 ).asInt();
	   double WaitTime = simulations[index].get("WaitTime", 0.1 ).asDouble();
	   double stim	 =	simulations[index].get("stimulus",   0 ).asDouble();
	   double s1dur	 =	simulations[index].get("s1duration", 0.10 ).asDouble();
	   int ncellPace =	simulations[index].get("ncellspaced", 1 ).asInt();
	   string blocker = simulations[index].get("blocker", "").asString();
	   double drug	 =	simulations[index].get("drug", 0e-6 ).asDouble();
	   int counterUpd = simulations[index].get("counterUpdateFreq", 10 ).asInt();
	   int ptupd	 =	simulations[index].get("printtimeFreq", 10 ).asInt();
	   int pCycupd	 =	simulations[index].get("printCycleFreq", 10 ).asInt();
	   int pMStAfter =	simulations[index].get("printMarkovStateAfter", MaxCycles-2  ).asInt();
	   bool logMark	 =	simulations[index].get("logMarkovIterations", false ).asBool();

	   s->set_params( dt, BCL, MaxCycles, WaitTime, stim, s1dur,
			   ncellPace, blocker, drug, counterUpd, ptupd, pCycupd,
			   pMStAfter, logMark);
	}
}


/*
 * readInputParameters - read input parameters from json file into sim list
 */
void readInputParameters( string inputfile, std::vector<simulation*> sims ){
	// Create an input file stream object to read in the
	//  text of the inputfile and then parse as a json file
	std::ifstream ifs(inputfile);
	string document((std::istreambuf_iterator<char>(ifs)),
						  std::istreambuf_iterator<char>());
	Json::Reader reader;
	Json::Value  root;
	bool parsedSuccessful = reader.parse( document, root );
	if (!parsedSuccessful){
		cerr << "Not able to parse input file!" << endl;
		return;
	}

	// At this point, the input has been parsed, we can start
	//  extracting all the data.
	const int nsims = root.get("NumberOfSimulations" , "0" ).asInt();
	if (nsims <= 0)
		return;

	std::string folder = root.get("OutputFolder" , "" ).asString();
	if (folder.compare( "" ) == 0)
		folder = "."; // Output will go into the current folder - project dir
	else {
		// Alternative to stuff below:
		// mkdir( str) - POSIX
		// or _mkdir( str ) in Win32 -- #include direct.h, stdlib.h
#if defined(_WIN32)
		string command = "mkdir " + folder;
#else
		string command = "mkdir -p " + folder;
#endif
		system(  command.c_str() );
	}
	cout << "Read input file; "  << nsims
			<< " simulations, output will go into: "<< folder << endl;

	const Json::Value simulations = root["Simulations"];
	const Json::Value cellparams  = root["Cells"];
	// Set cell parameters first - need a better way to set parameters.
	double params[3];
	params[0] = cellparams.get("GNa", "1.0").asDouble();
	params[1] = cellparams.get("GK1", "1.0").asDouble();
	params[2] = cellparams.get("Diffusion", "1.0").asDouble();
	for ( unsigned int index = 0; index < nsims; ++index ){
		// Iterates over the sequence elements.

	   int sim_type = simulations[index].get("type", "0" ).asInt();
	   int ncells = simulations[index].get("ncells", "1" ).asInt();
	   string fprefix = simulations[index].get("outputfileprefix", "").asString();

	   sims[index] = new simulation(sim_type, ncells, params, folder, fprefix);
	   simulation* s = sims[index];


		// Set up solver
	   string solvername = simulations[index].get("solver", "fe").asString();
	   solver* slvr = nullptr;
	   if (solvername.compare("fe") == 0){
			slvr = new solvers::forward_euler( s->mycells );
	   } else if (solvername.compare("rk2") == 0){
			slvr = new solvers::rk2( s->mycells );
	   } else if (solvername.compare("rk45") == 0){
			slvr = new solvers::rk45( s->mycells );
	   } else if (solvername.compare("rk45ad") == 0) {
		   slvr  = new solvers::rk54adaptive( s->mycells,
											1.0e-8,	// epsmin
											1.0e-4, // epsmax
											1.0e-6, // hmin
											1.0e-2  // hmax
		   );
	   } else if (solvername.compare("rl") == 0) {
		   slvr = new solvers::rushlarsen( s->mycells );
	   } else if (solvername.compare("rl2") == 0) {
		   slvr = new solvers::rushlarsen2( s->mycells );
	   }
	   s->setSolver( slvr );

	   double dt	 =	simulations[index].get("dt", 1e-6 ).asDouble();
	   double BCL	 =	simulations[index].get("BCL", 1000.0 ).asDouble();
	   int MaxCycles =	simulations[index].get("MaxCycles", 1 ).asInt();
	   double WaitTime = simulations[index].get("WaitTime", 0.1 ).asDouble();
	   double stim	 =	simulations[index].get("stimulus",   0 ).asDouble();
	   double s1dur	 =	simulations[index].get("s1duration", 0.10 ).asDouble();
	   int ncellPace =	simulations[index].get("ncellspaced", 1 ).asInt();
	   string blocker = simulations[index].get("blocker", "").asString();
	   double drug	 =	simulations[index].get("drug", 0e-6 ).asDouble();
	   int counterUpd = simulations[index].get("counterUpdateFreq", 10 ).asInt();
	   int ptupd	 =	simulations[index].get("printtimeFreq", 10 ).asInt();
	   int pCycupd	 =	simulations[index].get("printCycleFreq", 10 ).asInt();
	   int pMStAfter =	simulations[index].get("printMarkovStateAfter", MaxCycles-2  ).asInt();
	   bool logMark	 =	simulations[index].get("logMarkovIterations", false ).asBool();

	   s->set_params( dt, BCL, MaxCycles, WaitTime, stim, s1dur,
			   ncellPace, blocker, drug, counterUpd, ptupd, pCycupd,
			   pMStAfter, logMark);
	}
}
