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

	// set population manager pointer to nullptr immediately
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

    if ( parser -> GetDebuggerFlag() )
        Debug();
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
	bool analysis  = parser -> GetAnalysisFlag();
    std::size_t N  = parser -> GetNumParticles();

	// greet the user
	if (verbose) std::cout
        << "\n Building " << trials
        << " population(s) of size " << N << " ...\n";

    // iterate over all trials
    for (int t = 0; t < trials; t++){

        // display progress bar
        if (verbose == 2)
            display -> Progress(t, trials);

        // build a new population
        population -> Build(t);

        if (analysis){
//
            // find the nearest neighbor seperations
            population -> FindNeighbors(t);

            // fit a profile to curve
            population -> ProfileFit(t);
        }
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

void Simulation::Debug(){

    // FIXME: clean up notation for debugging output!!!!!
    // --------------------------------------------------

    // items for debugging mode
    std::vector<double> xlim = parser -> GetXlimits();
    std::vector<double> ylim = parser -> GetYlimits();
    std::vector<double> zlim = parser -> GetZlimits();

    std::vector<std::string> axes = parser -> GetAxes();
    std::vector<std::size_t> res  = parser -> GetResolution();

    std::string analysis_string = "";
    for ( int i = 0; i < axes.size(); i++ ){

        std::stringstream convert;

        if ( !analysis_string.empty() )
            convert << ", ";

        convert << "`" << axes[i] << "` (" << res[i] << ")";

        analysis_string += convert.str();
    }

    if ( !(parser -> GetAnalysisFlag()) )
        analysis_string = "None";

    std::string keep_pos    = parser -> GetKeepPosFlag()  ? "true" : "false";
    std::string keep_raw    = parser -> GetKeepRawFlag()  ? "true" : "false";
    std::string do_analysis = parser -> GetAnalysisFlag() ? "true" : "false";

    std::stringstream buffer;

    buffer << parser -> GetMeanBandwidth();
    std::string m_bandwidth = buffer.str();
    if ( m_bandwidth == "0" )
        m_bandwidth = "None";

    buffer.clear();
    buffer << parser -> GetStdevBandwidth();
    std::string s_bandwidth = buffer.str();
    if ( s_bandwidth == "0" )
        s_bandwidth = "None";

    std::cout << "\n" <<
    " Debugging Mode (on) | The following parameters are in use: \n"
    " -----------------------------------------------------------\n"
    "\n Number of Particles    = " << parser -> GetNumParticles() <<
    "\n Number of Trials       = " << parser -> GetNumTrials()    <<
    "\n Number of Threads      = " << parser -> GetNumThreads()   <<
    "\n Verbosity              = " << parser -> GetVerbosity()    <<
    "\n Keep Positions         = " << keep_pos <<
    "\n Keep Raw files         = " << keep_raw <<
    "\n X-limits               = (" << xlim[0] << ", " << xlim[1] << ")" <<
    "\n Y-limits               = (" << ylim[0] << ", " << ylim[1] << ")" <<
    "\n Z-limits               = (" << zlim[0] << ", " << zlim[1] << ")" <<
    "\n First Seed             = " << parser -> GetFirstSeed() <<
    "\n Sample Rate            = " << parser -> GetSampleRate() <<
    "\n Mean Bandwidth         = " << m_bandwidth <<
    "\n Stdev Bandwidth        = " << s_bandwidth <<
    "\n Analysis               = " << analysis_string <<
    "\n Output file pattern    = " << parser -> GetOutPath() << "*.dat" <<
    "\n Raw file pattern       = " << parser -> GetRawPath() << "*.dat" <<
    "\n Temporary file pattern = " << parser -> GetTmpPath() << "*.dat" <<
    "\n Position file pattern  = " << parser -> GetPosPath() << "*.dat" <<
    "\n RC file used           = " << parser -> GetRCFile() <<
    "\n" <<
    "\n Used PDFs:" <<
    "\n\n";

    for ( const auto& pdf : parser -> GetUsedPDFs() ){

        std::string pdftype = pdf.second.empty() ? "(Analytical)" :
            "(from file `" + pdf.second + "`)";

        std::cout << "\t * " << pdf.first << ", " << pdftype << std::endl;
    }
}

} // namespace Gaia
