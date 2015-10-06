// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/Monitor.hpp
//
// Header file for the Monitor class.
// This class contains tools and functions for displaying progress
// statements and `Monitor`ing computation time.
//

#ifndef _MONITOR_HH_
#define _MONITOR_HH_

#include <string>

namespace Gaia {

class Monitor {

public:

	static Monitor* GetInstance( );
	static void Release( );
	~Monitor( ) { }

	// initialization
	void Start( double freq = 0.25, double bar_width = 35.0 );

	// display progress bar
	void Progress( long iter, long max_iter, int threads = 1, 
		long first_iter = 0 );
	
	// display total time elapsed
	void ElapsedTime( );
	void TotalElapsedTime( );

	// construct local time of completion statement
	std::string CompletionTime( int );

private:

	static Monitor* instance;
	static bool started;
	Monitor( ) { }

	// timing and progress bar variables
	double _start, _last, _elapsed, _local_start_time;
	double _freq, _bar_width;

};

} // namespace Gaia 

#endif