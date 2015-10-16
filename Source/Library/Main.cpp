// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Library/Main.cc
//
// This is the `main` source file for the project. It simply makes the call
// to create the `Simulation` object and runs the code.

#include <iostream>
#include <exception>

#include "../Include/Simulation.hpp"
#include "../Include/Exception.hpp"

int main( const int argc, const char *argv[] ){

	try {

		// create and run the simulation
		Gaia::Simulation simulation(argc, argv);
		simulation.Run();

		return 0;

	} catch (const Gaia::Usage& usage){

		std::cout << usage.what() << std::endl;
		return 0;

	} catch (const Gaia::Exception& error){

		std::cerr << error.what() << std::endl;
		return 1;

	} catch (const std::exception& error){

		std::cerr << error.what() << std::endl;
		return 2;
	}
}
