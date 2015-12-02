// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// GNU General Public License v3.0
// Library/PopulationManager.cc
//
// #TODO:20 source

#include <omp.h>

#include <fstream>
#include <cmath>

#include <PopulationManager.hpp>
#include <ProfileManager.hpp>
#include <FileManager.hpp>
#include <Exception.hpp>
#include <KernelFit.hpp>
#include <Vector.hpp>
#include <Parser.hpp>
#include <Random.hpp>

#ifndef pi
#define pi 3.141592653589793
#endif

namespace Gaia {

PopulationManager::PopulationManager(){

	// initialize pointers to nullptr
	profiles  = nullptr;
	generator = nullptr;
}

PopulationManager::~PopulationManager()
{
	// delete profile manager
	if (profiles)
	{
		delete profiles;
		profiles = nullptr;
	}

	// delete PRNG
	if (generator)
	{
		delete generator;
		generator = nullptr;
	}
}

// set up the `Profile`s
void PopulationManager::Initialize(){

	// initialize profile manager
	profiles = new ProfileManager();
	profiles -> Initialize();

	// grab the file manager
	file = FileManager::GetInstance();

	// grab the Monitor
	display = Monitor::GetInstance();

	// grab the parser
	parser = Parser::GetInstance();

	// read in simulation parameters
	N          = parser -> GetNumParticles();
	first_seed = parser -> GetFirstSeed();
	threads    = parser -> GetNumThreads();
	trials     = parser -> GetNumTrials();
	verbose    = parser -> GetVerbosity();
	analysis   = parser -> GetAnalysisFlag();
	samples    = parser -> GetSampleRate() * double(N);

	mean_bandwidth  = parser -> GetMeanBandwidth();
	stdev_bandwidth = parser -> GetStdevBandwidth();

	// get cartesian limits for `box`
	Xlimits = parser -> GetXlimits();
	Ylimits = parser -> GetYlimits();
	Zlimits = parser -> GetZlimits();

	// initialize parallel mt19937 PRNG array
	generator = new ParallelMT(threads, first_seed);

	// initialize `positions` vector
	std::vector<Vector> new_population_vector;
	positions = new_population_vector;
	positions.resize(N);

	// build intervals on `positions` vector
	interval = Interval::Build(positions, threads);

	// compute the `span` of the space
	Vector span(
		Xlimits[1] - Xlimits[0], // maximum distance in `x`
		Ylimits[1] - Ylimits[0], // maximum distance in `y`
		Zlimits[1] - Zlimits[0]  // maximum distance in `z`
	);

	// initialization for FindNeighbors() algorithm
	max_seperation = span.Mag();

	// find furthest distance in each direction
	double max_X = std::abs(Xlimits[0]) < std::abs(Xlimits[1]) ?
		std::abs(Xlimits[1]) : std::abs(Xlimits[0]);

	double max_Y = std::abs(Ylimits[0]) < std::abs(Ylimits[1]) ?
		std::abs(Ylimits[1]) : std::abs(Ylimits[0]);

	double max_Z = std::abs(Zlimits[0]) < std::abs(Zlimits[1]) ?
		std::abs(Zlimits[1]) : std::abs(Zlimits[0]);

	double max_R   = std::sqrt(max_X*max_X + max_Y*max_Y);
	double max_Rho = std::sqrt(max_X*max_X + max_Y*max_Y + max_Z*max_Z);

	std::vector<double> Rlimits     = { 0.0, max_R    };
	std::vector<double> Rholimits   = { 0.0, max_Rho  };
	std::vector<double> Philimits   = { 0.0, 2.0 * pi };
	std::vector<double> Thetalimits = { 0.9, pi       };

	// build appropriate linespace limits for each axis
	std::map<std::string, std::vector<double>> Limits;
	Limits["X"    ] = Xlimits;
	Limits["Y"    ] = Ylimits;
	Limits["Z"    ] = Zlimits;
	Limits["R"    ] = Rlimits;
	Limits["Rho"  ] = Rholimits;
	Limits["Phi"  ] = Philimits;
	Limits["Theta"] = Thetalimits;

	// set the appropriate line-spaces for the analysis
	axis       = parser -> GetAxes();
	resolution = parser -> GetResolution();
	for ( int i = 0; i < axis.size(); i++ ){

		Axis[ axis[i] ] = Linespace(

			Limits[ axis[i] ][0], // start
			Limits[ axis[i] ][1], // end
			resolution[i]         // number of points for this axis
		);
	}

	// save map information to file
	file -> SaveMap( Axis );

	if (resolution.size() == 1){
		//FIXME: pooled mean/variance 1D initialization

		// initialize 1D pooled statistics vectors
		std::vector<double> init_1D(resolution[0], 0.0);

		pooled_mean_1D = init_1D;
		pooled_variance_1D = init_1D;

	} else if (resolution.size() == 2){
		//FIXME: pooled mean/variance 2D initialization

		// initialize 2D pooled statistics vectors
		std::vector< std::vector<double> > init_2D(resolution[0],
			std::vector<double>(resolution[1], 0.0));

		pooled_mean_2D = init_2D;
		pooled_variance_2D = init_2D;
	}

	// build coordinate function map for transposing `positions`
	Coord["X"] = X;
	Coord["Y"] = Y;
	Coord["Z"] = Z;
	Coord["R"] = R;
	Coord["Rho"] = Rho;
	Coord["Phi"] = Phi;
	Coord["Theta"] = Theta;
}

// build a new population set
void PopulationManager::Build(const int trial){

	if ( verbose > 2 ) std::cout
		<< "\n --------------------------------------------------"
		<< "\n Building population #" << trial + 1 << std::endl;

	#pragma omp parallel for
	for (int i = 0; i < threads; i++)
	for (std::size_t j = interval[i].start; j <= interval[i].end; j++){

        if ( verbose > 2 && !omp_get_thread_num() )
            display -> Progress(j, N, omp_get_num_threads() );

		// keep generating positions until we are `successful`
		while ( true ){

			// flag for determining condition for a `break`
			bool successful = true;

			// the new position vector (uniform in the `box`)
			Vector new_position(
				generator -> RandomReal( i, Xlimits ),
				generator -> RandomReal( i, Ylimits ),
				generator -> RandomReal( i, Zlimits ));

			// loop through PDFs and reject if less than uniform random number
			for ( const auto& pdf : profiles -> UsedPDFs ){

				ProfileBase *this_pdf = pdf;

				if ( this_pdf -> Evaluate(new_position) <
					generator -> RandomReal(i) ){

					successful = false;
					break;
				}
			}

			if (successful) {

				// keep the new position vector
				positions[j] = new_position;
				break;
			}
		}
	}

    if (verbose > 2)
        display -> Progress(N, N);

	// save results
	if ( parser -> GetKeepPosFlag() )
		file -> SavePositions(positions, trial + 1);
}

// solve for the nearest neighbor seperations
void PopulationManager::FindNeighbors(const int trial){

    if ( verbose )
        std::cout << "\n Computing seperations ..." << std::endl;
        std::cout.flush();

    std::vector<double> init(samples, max_seperation);
    seperations = init;

    #pragma omp parallel for
    for (std::size_t i = 0; i < samples; i++) {

        if ( verbose > 2 && !omp_get_thread_num() )
            display -> Progress(i, samples, omp_get_num_threads() );

        for (std::size_t j = 0; j < N; j++)
        if ( i != j ){

            double r = (positions[i] - positions[j]).Mag();

            if ( r < seperations[i] )
                seperations[i] = r;
        }
    }

    if ( verbose > 2 )
        display -> Progress(N, N);

    // save results
    if ( parser -> GetKeepRawFlag() )
        file -> SaveRaw(seperations, trial + 1);
}

// fit a curve/surface to the data from FindNeighbors()
void PopulationManager::ProfileFit(const int trial){

    // switch for 1D or 2D analysis, build KernelFit objects
    if ( Axis.size() == 1 ){

        if (verbose) std::cout
            << "\n Transposing vectors ... ";
            std::cout.flush();

        // build vector of coordinates (chosen at runtime)
        std::vector<double> coords(samples, 0.0);
        for ( std::size_t i = 0; i < samples; i++ )
            coords[i] = Coord[ axis[0] ]( positions[i] );

        if (verbose) std::cout
            << "done\n Solving profile with KernelFit1D ... \n";
            std::cout.flush();

        // initialize the KernelFit object
        KernelFit1D<double> kernel(coords, seperations, mean_bandwidth);

        // solve for the profile through the data
        std::vector<double> mean = kernel.Solve( Axis[ axis[0] ] );

        // set new bandwidth
        kernel.SetBandwidth(stdev_bandwidth);

        if (verbose)
		std::cout << "\n Solving for sample variances ... \n";
		std::cout.flush();

        // solve for the standard deviation of the fit
        std::vector<double> variance = kernel.Variance(Axis[ axis[0] ], true);

		// add results to cummulative results
		for ( std::size_t i = 0; i < resolution[0]; i++ ){

			pooled_mean_1D[i]     += mean[i];
			pooled_variance_1D[i] += variance[i];
		}

        // save the results to a file
        file -> SaveOutput(mean, variance, trial + 1);

    } else if ( Axis.size() == 2 ) {

		if (verbose) std::cout
			<< "\n Transposing vectors ... ";
			std::cout.flush();

		// build vector of coordinates (chosen at runtime)
		std::vector<double> coords_1(samples, 0.0);
		std::vector<double> coords_2(samples, 0.0);
		for ( std::size_t i = 0; i < samples; i++ ){

			coords_1[i] = Coord[ axis[0] ]( positions[i] );
			coords_2[i] = Coord[ axis[1] ]( positions[i] );
		}

		if (verbose) std::cout
			<< "done\n Solving profile with KernelFit2D ... \n";
			std::cout.flush();

		// initialize the KernelFit object
		KernelFit2D<double> kernel(coords_1, coords_2, seperations,
			mean_bandwidth);

		// solve for the profile through the data
		std::vector< std::vector<double> > mean = kernel.Solve(Axis[ axis[0] ],
			Axis[ axis[1] ] );

		// set new bandwidth
		kernel.SetBandwidth(stdev_bandwidth);

		if (verbose < 3)
			std::cout << "done";

		if (verbose)
			std::cout << "\n Solving for sample variances ... ";
			std::cout.flush();

		// solve for the variance of the fit
		std::vector< std::vector<double> > variance = kernel.Variance(
			Axis[ axis[0] ], Axis[ axis[1] ], true);

		if (verbose < 3)
			std::cout << "done";

		// add results to cummulative results
		for ( std::size_t i = 0; i < resolution[0]; i++ )
		for ( std::size_t j = 0; j < resolution[1]; j++ ){

			pooled_mean_2D[i][j]     += mean[i][j];
			pooled_variance_2D[i][j] += variance[i][j];
		}

		// save the results to a file
		file -> SaveOutput(mean, variance, trial + 1);

    } else throw Exception("\n Error: From PopulationManager::ProfileFit, "
        "something is wrong. Axis.size() > 2");
}

// combine statistics for all trials
void PopulationManager::Analysis(){

	// take the mean of the pooled statistics and output results

	if (verbose)
		std::cout << "\n Pooling statistics ... ";
		std::cout.flush();

	// switch for 1D or 2D analysis, build KernelFit objects
	if ( Axis.size() == 1 ){

		for (auto &x : pooled_mean_1D)
			x /= trials;

		for (auto &x : pooled_variance_1D)
			x /= trials;

		if (verbose)
			std::cout << "done";

		file -> SaveOutput(pooled_mean_1D, pooled_variance_1D, 0);

    } else if ( Axis.size() == 2 ) {

		for (std::size_t i = 0; i < resolution[0]; i++)
		for (std::size_t j = 0; j < resolution[1]; j++){

			pooled_mean_2D[i][j]     /= trials;
			pooled_variance_2D[i][j] /= trials;
		}

		if (verbose)
			std::cout << "done";

		file -> SaveOutput(pooled_mean_2D, pooled_variance_2D, 0);

	} else throw Exception("\n Error: From PopulationManager::Analysis, "
        "something is wrong. Axis.size() > 2");
}

std::vector<double> PopulationManager::Linespace(const double start,
    const double end, const std::size_t length){

    // initialize the vector
    std::vector<double> linespace( length, 0.0 );

    linespace[0] = start;
    double dx    = (end - start) / double(length - 1);

    for ( std::size_t i = 1; i < length; i++ )
        linespace[i] = linespace[i-1] + dx;

    return linespace;
}

std::vector<Interval> Interval::Build(const std::vector<Vector> &input,
	const std::size_t num){
	//
	// `Build` a vector of `Interval`s based on an `input` vector and the
	// `num`ber of subdivisions requested.
	//

	// initialize the output vector
	std::vector<Interval> output;

	// the size of intervals
	std::size_t interval_size = input.size() / num;

	// the first interval starts with the first element of the `input`
	// the name `previous` will make sense further down
	Interval previous(0, interval_size);
	output.push_back(previous);

	// leap-frog through the intervals
	for (std::size_t i = 0; i < num - 1; i++){

		Interval next(previous.end + 1, previous.end + interval_size);
		output.push_back(next);

		previous = next;
	}

	// rectify final interval (for odd lengths)
	output[num-1].end = input.size() - 1;

	return output;
}

} // namespace Gaia
