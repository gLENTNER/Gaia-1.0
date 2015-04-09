// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Library/PopulationManager.cc
//
// TODO: source

#include <PopulationManager.hh>

namespace GAIA {

PopulationManager::PopulationManager(){
	
	// initialize pointers to NULL
	profiles = generator = positions = NULL;
	
	// initialize profile manager
	profiles = new ProfileManager();
	
	// grab the file manager
	file = FileManager::GetInstance();
	
	// read in simulation parameters
	unsigned long long N          = parser -> GetNumParticles();
	unsigned long long first_seed = parser -> GetFirstSeed(); 
	int                threads    = parser -> GetNumThreads();
	
	// get cartesian limits for `box`
	std::vector<double> _xlimits = parser -> GetXlimits();
	std::vector<double> _ylimits = parser -> GetYlimits();
	std::vector<double> _zlimits = parser -> GetZlimits();
	
	// initialize parallel mt19937 PRNG array
	generator = new ParallelMT(threads, first_seed);
	
	// create vector of `Vector`s for positions
	positions = new std::vector<Vector>;
	positions -> resize(N);
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
	
	// delete positions
	if (positions){
		delete positions;
		positions = NULL;
	}
}

// set up the `Profile`s
void PopulationManager::Initialize(){
	
	// initialize the `Profile`s in the profile manager
	profiles -> Initialize();
}

// build a new population set
void PopulationManager::Build(const int trial){
	
	#pragma omp parallel for
	for (int i = 0; i < threads; i++)
	for (unsigned long long j = interval[i].start; j <= interval[i].end; j++){
		
		// keep generating positions until we are `successful`
		while ( true ){
			
			// flag for determining condition for a `break`
			bool successful = true;
			
			// the new position vector (uniform in the `box`)
			Vector new_position(
				generator -> RandomReal( i, _xlimits ),
				generator -> RandomReal( i, _ylimits ),
				generator -> RandomReal( i, _zlimits )
			);
			
			// loop through PDFs and reject if less than uniform random number
			for (std::list<ProfileBase*> pdf = profiles -> usedPDFs.begin();
				pdf != profiles -> usedPDFs.end(); pdf++){
			
				ProfileBase *this_pdf = *pdf;
				
				if ( this_pdf -> Evaluate(new_position) < 
					generator -> RandomReal() ){
					
					successful = false; break;	
				}
				
			} if (successful) {
				// keep the new position vector
				positions[j] = new_position;
				break;
			}
		}
	}
	
	// save results
	if ( parser -> GetKeepPosFlag() )
		file -> SavePositions(positions, trial);
}

// solve for the nearest neighbor separations
void PopulationManager::FindNeighbors(const int trial){
	
}

// fit a curve/surface to the data from FindNeighbors()
void PopulationManager::ProfileFit(const int trial){
	
}

// combine statistics for all trials
void PopulationManager::Analysis(){
	
}


} // namespace GAIA