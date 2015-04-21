// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Library/Simulation.cc
//
// This source file gives the definitions for the `Simulation` class. This
// object is the manager for the GAIA simulation. It essentially just puts all
// the pieces of the code together.

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <Simulation.hh>
#include <FileManager.hh>
#include <Monitor.hh>
#include <PopulationManager.hh>
#include <Parser.hh>

namespace Gaia {

Simulation::Simulation(const int argc, const char *argv[]){

	// set population manager pointer to NULL immediately
	population = nullptr;

	// start the clock right away
	display = Monitor::GetInstance();
	display -> Start();

	// create parser and interpret command line arguments + RC file
	parser = Parser::GetInstance();
	parser -> Setup(argc, argv);
    
    if ( parser -> GetVerbosity() ) std::cout
        << "\n Welcome to GAIA | version 1.0.1"
        << "\n Copyright (c) Geoffrey Lentner 2015 (GPLv3)\n"
        << "\n";

	// create and initialize the population manager
    population = new PopulationManager();
    population -> Initialize();

	// // create and initialize the file manager
    file = FileManager::GetInstance();
    file -> Initialize();

	// FIXME: clean up notation for debugging output!!!!!
	// --------------------------------------------------

	// items for debugging mode
	std::vector<double> xlim = parser -> GetXlimits();
	std::vector<double> ylim = parser -> GetYlimits();
	std::vector<double> zlim = parser -> GetZlimits();
	std::string xyres;
	int radial_res;
	std::string analysis_coord;
	int analysis = parser -> GetAnalysisDomain();
	std::stringstream convert;
	std::vector<int> xy;
	switch ( analysis ){
		case 1:
			analysis_coord  = "R (Cyclindrical)";
			radial_res = parser -> GetRadialResolution();
			break;
		case 2:
			analysis_coord  = "Rho (Spherical)";
			radial_res = parser -> GetRadialResolution();
			break;
		case 3:
			analysis_coord  = "XY (Plane)";
			xy = parser -> GetXYResolution();
			convert << "(" << xy[0] << ", " << xy[1] << ")\n";
			xyres = convert.str();
			break;

		default:
			analysis_coord = "None (given --no-analysis)";
	}

	std::string keep_pos    = parser -> GetKeepPosFlag()  ? "true" : "false";
	std::string keep_raw    = parser -> GetKeepRawFlag()  ? "true" : "false";
	std::string do_analysis = parser -> GetAnalysisFlag() ? "false" : "true";

	if (analysis != 3) xyres = "None (unspecified)";
	convert.clear();
	std::string rres;
	if ( analysis == 1 || analysis == 2 ) {
		convert << radial_res;
		rres = convert.str();

	} else rres = "None (unspecified)";

	// show parameters if in `debug` mode
	if ( parser -> GetDebuggerFlag() ){
		std::cout << "\n" <<
		" Debugging Mode (on) | The following parameters are in use: \n"
		" -----------------------------------------------------------\n"
		"\n Number of Particles    = " << parser -> GetNumParticles() <<
		"\n Number of Trials       = " << parser -> GetNumTrials()    <<
		"\n Number of Threads      = " << parser -> GetNumThreads()   <<
		"\n Verbosity              = " << parser -> GetVerbosity()    <<
		"\n Keep Positions         = " << keep_pos <<
		"\n Keep Raw files         = " << keep_raw <<
		"\n Do Analysis            = " << do_analysis <<
		"\n X-limits               = (" << xlim[0] << ", " << xlim[1] << ")" <<
		"\n Y-limits               = (" << ylim[0] << ", " << ylim[1] << ")" <<
		"\n Z-limits               = (" << zlim[0] << ", " << zlim[1] << ")" <<
		"\n Analysis Domain        = " << analysis_coord <<
		"\n X-Y Resolution         = " << xyres <<
		"\n Radial Resolution      = " << radial_res <<
		"\n Output file pattern    = " << parser -> GetOutPath() << "*.dat" <<
		"\n Raw file pattern       = " << parser -> GetRawPath() << "*.dat" <<
		"\n Temporary file pattern = " << parser -> GetTmpPath() << "*.bin" <<
		"\n Position file pattern  = " << parser -> GetPosPath() << "*.dat" <<
		"\n RC file used           = " << parser -> GetRCFile() <<
		"\n First Seed             = " << parser -> GetFirstSeed() <<
		"\n" <<
		"\n Used PDFs:" <<
		"\n -----------------------\n\n";

		for ( const auto& pdf : parser -> GetUsedPDFs() ){
			std::cout << " " << pdf.first;
			std::cout << " `" << pdf.second << "`\n";
		}
	}
}

Simulation::~Simulation(){

	// release the parser
	Parser::Release();

	// release file manager
	FileManager::Release();

	// release the display monitor
	Monitor::Release();

	// delete population manager
    if (population){
        delete population;
        population = nullptr;
    }
}

void Simulation::Run(){

	// get simulation parameters
	int verbose    = parser -> GetVerbosity();
	int trials     = parser -> GetNumTrials();
	int analysis   = parser -> GetAnalysisDomain();
    std::size_t  N = parser -> GetNumParticles();

	// greet the user
	if (verbose) std::cout << "\n Building " << trials <<
        " population(s) of size " << N << " ...\n";

    // iterate over all trials
    for (int t = 0; t < trials; t++){

        // display progress bar
        if (verbose == 2)
            display -> Progress(t, trials);

        // build a new population
        population -> Build(t);

//        if (analysis){
//
//            // find the nearest neighbor seperations
//            population -> FindNeighbors(t);
//
//            // fit a profile to curve
//            population -> ProfileFit(t);
//        }
    }

    // complete progress bar
    if (verbose == 2)
        display -> Progress(trials, trials);

    // combine statistics for nearest neighbor analysis
//    if (analysis)
//        population -> Analysis();

	if (verbose)
		display -> TotalElapsedTime();
}

} // namespace Gaia
