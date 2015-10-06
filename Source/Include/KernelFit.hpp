// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/KernelFit.hpp

// This header file contains the template for the KernelFit objects.


#ifndef _KERNELFIT_HH_
#define _KERNELFIT_HH_

#include <string>
#include <vector>
#include <cmath>

#include <Exception.hpp>
#include <Parser.hpp>
#include <Monitor.hpp>

namespace Gaia {

template<class T>
class KernelFit1D {

public:

	KernelFit1D(){}
	KernelFit1D(const std::vector<T> &x, const std::vector<T> &y,
		const T &bandwidth);

	// kernel function used by default
	T Kernel(const T &x){return exp( -0.5 * x * x / _b );}

	// solve for smooth curve through data
	std::vector<T> Solve(const std::vector<T> &x, const bool unbiased = false);

	// solve by alternative kernel function
	std::vector<T> Solve(const std::vector<T> &x, T (*W)(T),
        const bool unbiased = false);

	// solve for estimated deviations
    std::vector<T> Variance(const std::vector<T> &x,
		const bool unbiased = false);

    // solve for estimated deviations by alternative kernel function
    std::vector<T> Variance(const std::vector<T> &x, T (*W)(T),
        const bool unbiased = false);

	// solve for estimated deviations
    std::vector<T> StdDev(const std::vector<T> &x, const bool unbiased = false);

    // solve for estimated deviations by alternative kernel function
    std::vector<T> StdDev(const std::vector<T> &x, T (*W)(T),
        const bool unbiased = false);

	// set with function to ensure it is squared
	void SetBandwidth(T bandwidth){ _b = bandwidth*bandwidth; }

protected:

	T _b;
    std::vector<T> _x, _y;
    std::size_t N;

	Parser *parser;
	Monitor *display;
	int verbose;
};

template<class T>
class KernelFit2D {

public:

	KernelFit2D(){}
	KernelFit2D(const std::vector<T> &x, const std::vector<T> &y,
		const std::vector<T> &z, const T &bandwidth);

	// kernel function used by default
	T Kernel(const T &x, const T &y){ return exp( -0.5 * (x*x + y*y) / _b); }

	// solve for the smooth surface through the data
	std::vector< std::vector<T> > Solve(const std::vector<T> &x,
		const std::vector<T> &y, const bool unbiased = false);

    // solve by alternative kernel function
    std::vector< std::vector<T> > Solve(const std::vector<T> &x,
        const std::vector<T> &y, T (*W)(T, T), const bool unbiased = false);

	// solve for estimated deviations
    std::vector< std::vector<T> > Variance(const std::vector<T> &x,
    	const std::vector<T> &y, const bool unbiased = false);

    // solve for estimated deviations
    std::vector< std::vector<T> > Variance(const std::vector<T> &x,
        const std::vector<T> &y, T (*W)(T, T), const bool unbiased = false);

    // solve for estimated deviations
    std::vector< std::vector<T> > StdDev(const std::vector<T> &x,
    	const std::vector<T> &y, const bool unbiased = false);

    // solve for estimated deviations
    std::vector< std::vector<T> > StdDev(const std::vector<T> &x,
        const std::vector<T> &y, T (*W)(T, T), const bool unbiased = false);


	// set with function to ensure it is squared
	void SetBandwidth(T bandwidth){ _b = bandwidth*bandwidth; }

protected:

	T _b;
	std::vector<T> _x, _y, _z;
    std::size_t N;

	Parser *parser;
	Monitor *display;
	int verbose;
};

// exception thrown by KernelFit objects
class KernelFitError : public Exception {
public:

	KernelFitError(const std::string& msg): Exception(
		"\n --> KernelFitError: " + msg){ }
};

} // namespace Gaia

#endif
