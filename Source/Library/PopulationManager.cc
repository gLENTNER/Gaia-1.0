// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Library/PopulationManager.cc
//
// TODO: source

#include <omp.h>

#include <fstream>
#include <cmath>

#include <PopulationManager.hh>
#include <ProfileManager.hh>
#include <FileManager.hh>
#include <Exception.hh>
#include <KernelFit.hh>
#include <Vector.hh>
#include <Parser.hh>
#include <Random.hh>

#ifndef pi
#define pi 3.141592653589793
#endif

namespace Gaia {

PopulationManager::PopulationManager(){

	// initialize pointers to nullptr
	profiles  = nullptr;
	generator = nullptr;
}

PopulationManager::~PopulationManager(){

	// delete profile manager
	if (profiles){
		delete profiles;
		profiles = nullptr;
	}

	// delete PRNG
	if (generator){
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

					successful = false; break;
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
            << "done\n Solving 1D KernelFit ... ";
            std::cout.flush();

        // initialize the KernelFit object
        KernelFit1D<double> kernel(coords, seperations, mean_bandwidth);

        // solve for the profile through the data
        std::vector<double> mean = kernel.Solve( Axis[ axis[0] ] );

        // set new bandwidth
        kernel.SetBandwidth(stdev_bandwidth);

        if (verbose) std::cout
            << "done\n Solving for standard deviations ... ";
            std::cout.flush();

        // solve for the standard deviation of the fit
        std::vector<double> stdev = kernel.StdDev( Axis[ axis[0] ] );

        if (verbose) std::cout
            << "done\n";
            std::cout.flush();

        // save the results to a file
        file -> SaveTemp(Axis[ axis[0] ], mean, stdev, trial + 1);

//    } else if ( Axis.size() == 2 ) {


    } else throw Exception("\n Error: From PopulationManager::ProfileFit, "
        "something is wrong. Axis.size() > 2");
}

// combine statistics for all trials
void PopulationManager::Analysis(){

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
	// the name `last` will make sense further down
	Interval last(0, interval_size);
	output.push_back(last);

	// leap-frog through the intervals
	for (std::size_t i = 0; i < num - 1; i++){

		Interval next(last.end + 1, last.end + interval_size);
		output.push_back(next);

		last = next;
	}

	// rectify final interval (for odd lengths)
	output[num-1].end = input.size() - 1;

	return output;
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

} // namespace Gaia
