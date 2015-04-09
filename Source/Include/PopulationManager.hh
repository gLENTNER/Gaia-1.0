// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/PopulationManager.hh
//
// TODO: header

#ifndef _POPULATIONMANAGER_HH_
#define _POPULATIONMANAGER_HH_

#include <FileManager.hh>
#include <ProfileManager.hh>
#include <Random.hh>
#include <Vector.hh>

namespace GAIA {

class PopulationManager {

public:
	
	PopulationManager();
	~PopulationManager();
	
	// set up the `Profile`s
	Initialize();
	
	// build a new population set
	Build(const int trial);
	
	// solve for the nearest neighbor separations
	FindNeighbors(const int trial);
	
	// fit a curve/surface to the data from FindNeighbors()
	ProfileFit(const int trial);
	
	// combine statistics for all trials
	Analysis(); 
	
private:
	
	// file manager is a private member
	FileManager *file;
	
	// profile manager keeps `Profile`s
	ProfileManager *profiles;
	
	// parallel mt19937 PRNG array
	ParallelMT *generator;
	
	// vector of `Vector` positions
	std::vector<Vector> *positions;

};

} // namespace GAIA

#endif