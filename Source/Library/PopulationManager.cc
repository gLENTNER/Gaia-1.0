// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Library/PopulationManager.cc
//
// TODO: source

#include <omp.h>

#include <fstream>

#include <PopulationManager.hh>
#include <ProfileManager.hh>
#include <FileManager.hh>
#include <Vector.hh>
#include <Parser.hh>
#include <Random.hh>

namespace Gaia {

PopulationManager::PopulationManager(){

	// initialize pointers to NULL
	profiles  = nullptr;
	generator = nullptr;
}

PopulationManager::~PopulationManager(){

	// delete profile manager
	if (profiles){
		delete profiles;
		profiles = NULL;
	}

	// delete PRNG
	if (generator){
		delete generator;
		generator = NULL;
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
    maximum_seperation = span.Mag();
}

// build a new population set
void PopulationManager::Build(const int trial){

    if ( verbose > 2 )
        std::cout << "\n Building population #" << trial << std::endl;
    
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
		file -> SavePositions(positions, trial);
}

// solve for the nearest neighbor seperations
void PopulationManager::FindNeighbors(const int trial){

    if ( verbose )
        std::cout << "\n Computing seperations ...";
    
    std::vector<double> init(N, maximum_seperation);
    seperations = init;
    
    #pragma omp parallel for
    for (std::size_t i = 0; i < N; i++)
    for (std::size_t j = 0; j < N; j++)
    if ( i != j ){
        
        double r = (positions[i] - positions[j]).Mag();
        
        if ( r < seperations[i] )
            seperations[i] = r;
    }
    
    if ( verbose )
        std::cout << " done";
    
    // save results
    if ( parser -> GetKeepRawFlag() )
        file -> SaveRaw(seperations, trial);
}

// fit a curve/surface to the data from FindNeighbors()
void PopulationManager::ProfileFit(const int trial){

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

} // namespace Gaia
