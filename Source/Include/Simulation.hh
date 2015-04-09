// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/Simulation.hh
//
// This header file gives the template for the `Simulation` class. This object
// is the manager for the GAIA simulation. It essentially just puts all the
// pieces of the code together.


#ifndef _SIMULATION_HH_
#define _SIMULATION_HH_

#include <Parser.hh>

namespace GAIA {

class Simulation {

public:
	
	Simulation(const int argc, const char *argv[]);
	~Simulation();
	
	void Run();
	
private:
	
	Parser *parser;
	Monitor *display;
	FileManager *file;
	PopulationManager *population;

};

} // namespace GAIA

#endif