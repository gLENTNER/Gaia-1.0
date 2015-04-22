// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/PopulationManager.hh
//
// TODO: header

#ifndef _POPULATIONMANAGER_HH_
#define _POPULATIONMANAGER_HH_

#include <vector>

#include <FileManager.hh>
#include <ProfileManager.hh>
#include <Monitor.hh>
#include <Random.hh>
#include <Vector.hh>

namespace Gaia {

// helper class to manager segments of vectors
class Interval {
public:

	Interval(std::size_t n, std::size_t m): start(n), end(m){}

	// only two members
	std::size_t start, end;

	// factory function returns vector of intervals
	static std::vector<Interval> Build(const std::vector<Vector> &input,
		const std::size_t num);
};

class PopulationManager {

public:

	PopulationManager();
	~PopulationManager();

	// set up the `Profile`s
	void Initialize();

	// build a new population set
	void Build(const int trial);

	// solve for the nearest neighbor separations
	void FindNeighbors(const int trial);

	// fit a curve/surface to the data from FindNeighbors()
	void ProfileFit(const int trial);

	// combine statistics for all trials
	void Analysis();

private:

	// file manager is a private member
	FileManager *file;

	// profile manager keeps `Profile`s
	ProfileManager *profiles;

	// parallel mt19937 PRNG array
	ParallelMT *generator;
    
    // Monitor
    Monitor *display;
    
    // parser
    Parser *parser;

	// limits for volume
	std::vector<double> Xlimits, Ylimits, Zlimits;

	// vector of `Vector` positions
	std::vector<Vector>   positions;
	std::vector<Interval> interval;

	// nearest neighbor seperations
	std::vector<double> seperations;
    double max_seperation;

	// simulation parameters from parser
	std::size_t N;
	unsigned long long first_seed;
	int threads, verbose;

};


} // namespace Gaia

#endif
