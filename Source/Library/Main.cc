// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Library/Main.cc
//
// This is the `main` source file for the project. It simply makes the call
// to create the `Simulation` object and runs the code.

#include <iostream>
#include <exception>

#include <Simulation.hh>
#include <Exception.hh>

int main( const int argc, const char *argv[] ){
	
	try {

		// create and run the simulation
		GAIA::Simulation simulation(argc, argv);
		simulation.Run();
		
	} catch (const GAIA::Usage& usage){
		
		std::cout << usage.what() << std::endl;
		return 0;
	
	} catch (const GAIA::Exception& error){
		
		std::cerr << error.what() << std::endl;
		return 1;
		
	} catch (const std::exception& error){
		
		std::cerr << error.what() << std::endl;
		return 2;
	}
	
	return 0;
}