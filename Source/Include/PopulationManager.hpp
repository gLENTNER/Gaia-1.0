// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/PopulationManager.hpp
//
// TODO: header

#ifndef _POPULATIONMANAGER_HH_
#define _POPULATIONMANAGER_HH_

#include <vector>

#include <FileManager.hpp>
#include <ProfileManager.hpp>
#include <Monitor.hpp>
#include <Random.hpp>
#include <Vector.hpp>

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

    // helper function for building the `Axis` map
    std::vector<double> Linespace(const double, const double, const std::size_t);

	// limits for volume
	std::vector<double> Xlimits, Ylimits, Zlimits;

    // map analysis coordinates to line-spaces
    std::map<std::string, std::vector<double>> Axis;
    std::vector<std::string> axis;
    std::vector<std::size_t> resolution;

    // map of functions, (just like in ProfileBase)
    std::map< std::string, double (*)(const Vector&) > Coord;

    // functions in the above map (calls to vector coordinates)
    static double X(const Vector &vec)     { return vec.X();     }
    static double Y(const Vector &vec)     { return vec.Y();     }
    static double Z(const Vector &vec)     { return vec.Z();     }
    static double R(const Vector &vec)     { return vec.R();     }
    static double Rho(const Vector &vec)   { return vec.Rho();   }
    static double Phi(const Vector &vec)   { return vec.Phi();   }
    static double Theta(const Vector &vec) { return vec.Theta(); }

	// vector of `Vector` positions
	std::vector<Vector>   positions;
	std::vector<Interval> interval;

	// nearest neighbor seperations and match to specified coordinates
	std::vector<double> seperations, pooled_mean_1D, pooled_variance_1D;
	std::vector< std::vector<double> > pooled_mean_2D, pooled_variance_2D;
    double max_seperation;

	// simulation parameters from parser
	std::size_t N, samples;
	unsigned long long first_seed;
	int threads, trials, verbose;
    bool analysis;
    double mean_bandwidth, stdev_bandwidth;

};


} // namespace Gaia

#endif
