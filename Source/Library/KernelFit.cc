// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Library/KernelFit.cc

// This source file contains the definitions for the KernelFit objects.

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>

#include <KernelFit.hh>
#include <Parser.hh>
#include <Monitor.hh>

namespace Gaia {

template<class T>
KernelFit1D<T>::KernelFit1D(const std::vector<T> &x, const std::vector<T> &y,
	const T &bandwidth){

	// Constructor for the KernelFit1D object. Save the `x`, `y` data
	// and set an initial value for the `bandwidth`.

	if ( x.empty() || y.empty() )
		throw KernelFitError("From KernelFit1D::KernelFit1D(), "
		"one or both input vectors are empty!");

	if ( x.size() != y.size() )
		throw KernelFitError("From KernelFit1D::KernelFit1D(), input vectors "
		"must be equal in length!");

	if ( bandwidth <= 0.0 )
		throw KernelFitError("From KernelFit1D::KernelFit1D(), the bandwidth "
		"must be greater than zero!");

	_x = x;
	_y = y;
	_b = bandwidth * bandwidth; // squared ahead of time
	N  = x.size(); // they're all the same length...

	// include display monitor for Gaia project
	parser = Parser::GetInstance();
	verbose = parser -> GetVerbosity();
	display = Monitor::GetInstance();

}

template<class T>
std::vector<T> KernelFit1D<T>::Solve(const std::vector<T> &x,
	const bool unbiased){

	//
	// solve for the smooth profile through the data at all `x`
	//

	if ( x.empty() )
        throw KernelFitError("From KernelFit1D::Solve(), the input vector "
        "cannot be empty!");

	std::vector<T> f( x.size(), 0.0);

	// omp_set_num_threads() should be called prior to here!
	#pragma omp parallel for shared(f)
	for (std::size_t i = 0; i <  x.size(); i++){

		if ( verbose > 2 && !omp_get_thread_num() )
            display -> Progress(i, x.size(), omp_get_num_threads() );

		T sum = 0.0;

		for (std::size_t j = 0; j < N; j++){

			T W   = Kernel(_x[j] - x[i]);
			f[i] += W * _y[j];
			sum  += W;
		}

		if (unbiased){

			// adjust the sum over weights to `unbias` the result
			// only relavent when called from Stdev()!!!
			f[i] /= (1.0 - 1.0 / N) * sum;

		} else f[i] /= sum;
	}

	if (verbose > 2)
		display -> Progress(1, 1); // complete

	return f;
}

template<class T>
std::vector<T> KernelFit1D<T>::Solve(const std::vector<T> &x, T (*W)(T),
	const bool unbiased){

    //
    // solve for the smooth profile through the data at all `x`
    // using an alternative kernel function `W`
    //

    if ( x.empty() )
        throw KernelFitError("From KernelFit1D::Solve(), the input vector "
        "cannot be empty!");

    std::vector<T> f( x.size(), 0.0);

    // omp_set_num_threads() should be called prior to here!
    #pragma omp parallel for shared(f)
    for (std::size_t i = 0; i <  x.size(); i++){

		if ( verbose > 2 && !omp_get_thread_num() )
            display -> Progress(i, x.size(), omp_get_num_threads() );

        T sum = 0.0;

        for (std::size_t j = 0; j < N; j++){

            T WW  = W(_x[j] - x[i]);
            f[i] += WW * _y[j];
            sum  += WW;
        }

		if (unbiased){

			// adjust the sum over weights to `unbias` the result
			// only relavent when called from Stdev()!!!
			f[i] /= (1.0 - 1.0 / N) * sum;

		} else f[i] /= sum;
    }

	if (verbose > 2)
		display -> Progress(1, 1); // complete

    return f;
}

template<class T>
std::vector<T> KernelFit1D<T>::Variance(const std::vector<T> &x,
	const bool unbiased){

	//
    // Solve for the estimated standard deviation by evaluating
    // the profile *at* the raw data points.
    //

    if ( x.empty() )
        throw KernelFitError("From KernelFit1D::Variance(), the input vector "
        "cannot be empty!");

    // solve profile at data points
    std::vector<T> f = Solve( _x );

    // solve variance at data points
    std::vector<T> var(N, 0.0);
    for (std::size_t i = 0; i < N; i++)
        var[i] = std::pow(_y[i] - f[i], 2.0);

    // solve for smooth curve through variance points
    KernelFit1D<T> profile(_x, var, _b);

    return profile.Solve(x, unbiased);
}

template<class T>
std::vector<T> KernelFit1D<T>::Variance(const std::vector<T> &x,
	T (*W)(T), const bool unbiased){

	//
    // Solve for the estimated variance by evaluating
    // the profile *at* the raw data points. In this version,
    // I use an alternative kernel function given by the user.
	//

    if ( x.empty() )
        throw KernelFitError("From KernelFit1D::Variance(), the input vector "
        "cannot be empty!");

    // solve profile at data points
    std::vector<T> f = Solve( _x, W );

    // solve variance at data points
    std::vector<T> var(_x.size(), 0.0);
    for (std::size_t i = 0; i < _x.size(); i++)
        var[i] = std::pow(_y[i] - f[i], 2.0);

    // solve for smooth curve through variance points
    KernelFit1D<T> profile(_x, var, _b);

    return profile.Solve(x, W, unbiased);
}

template<class T>
std::vector<T> KernelFit1D<T>::StdDev(const std::vector<T> &x,
	const bool unbiased){

	//
    // Solve for the estimated standard deviation by evaluating
    // the profile *at* the raw data points.
    //

	if ( x.empty() )
        throw KernelFitError("From KernelFit1D::StdDev(), the input vector "
        "cannot be empty!");

	std::vector<T> stdev = Variance(x, unbiased);

    // take sqrt for standard deviation
    for ( auto& x : stdev )
        x = std::sqrt(x);

    return stdev;
}

template<class T>
std::vector<T> KernelFit1D<T>::StdDev(const std::vector<T> &x,
	T (*W)(T), const bool unbiased){

	//
    // Solve for the estimated standard deviation by evaluating
    // the profile *at* the raw data points. In this version,
    // I use an alternative kernel function given by the user.
	//

	if ( x.empty() )
        throw KernelFitError("From KernelFit1D::StdDev(), the input vector "
        "cannot be empty!");

	std::vector<T> stdev = Variance(x, W, unbiased);

    // take sqrt for standard deviation
    for ( auto& x : stdev )
        x = std::sqrt(x);

    return stdev;
}

template<class T>
KernelFit2D<T>::KernelFit2D(const std::vector<T> &x, const std::vector<T> &y,
	const std::vector<T> &z, const T &bandwidth){

	// Constructor for the KernelFit1D object. Save the `x`, `y` data
	// and set an initial value for the `bandwidth`.

	if ( x.empty() || y.empty() || z.empty() )
		throw KernelFitError("From KernelFit2D::KernelFit2D(), one or more "
			"input vectors were empty!");

	if ( x.size() != y.size() || x.size() != z.size() )
		throw KernelFitError("From KernelFit2D::KernelFit2D(), input vectors "
			"must be equal in length!");

	if ( bandwidth <= 0.0 )
		throw KernelFitError("From KernelFit2D::KernelFit2D(), the bandwidth "
			"must be greater than zero!");

	_x = x;
	_y = y;
	_z = z;
	_b = bandwidth * bandwidth; // square
	N  = x.size(); // they're all the same length...

	// include display monitor for Gaia project
	parser = Parser::GetInstance();
	verbose = parser -> GetVerbosity();
	display = Monitor::GetInstance();
}

template<class T>
std::vector< std::vector<T> > KernelFit2D<T>::Solve(const std::vector<T> &x,
	const std::vector<T> &y, const bool unbiased){

	//
	// solve for the smooth surface through the data at all (x, y)
	//

	if ( x.empty() || y.empty() )
		throw KernelFitError("From KernelFit2D::Solve(), one or both of "
			"`x` and `y` were empty!");

	// initialize f[x][y] to zeros with proper dimensions
	std::vector< std::vector<T> > f(x.size(), std::vector<T>(y.size(), 0.0));

	// omp_set_num_threads() should be called prior to here!
	#pragma omp parallel for shared(f)
	for (std::size_t i = 0; i < x.size(); i++){

		if ( verbose > 2 && !omp_get_thread_num() )
			display -> Progress(i, x.size(), omp_get_num_threads() );

		for (std::size_t j = 0; j < y.size(); j++){

			T sum = 0.0;

			for (std::size_t k = 0; k < N; k++){

				T W      = Kernel(x[i] - _x[k], y[j] - _y[k]);
				f[i][j] += W * _z[k];
				sum     += W;
			}

			if (unbiased){

				// adjust the sum over weights to `unbias` the result
				// only relavent when called from Stdev()!!!
				f[i][j] /= (1.0 - 1.0 / N) * sum;

			} else f[i][j] /= sum;
		}
	}

	if (verbose > 2)
		display -> Progress(1, 1); // complete

	return f;
}

template<class T>
std::vector< std::vector<T> > KernelFit2D<T>::Solve(const std::vector<T> &x,
    const std::vector<T> &y, T (*W)(T, T), const bool unbiased){

    //
    // solve for the smooth surface through the x,y data using alternative
    // kernel function `W`.
    //

    if ( x.empty() || y.empty() )
        throw KernelFitError("From KernelFit2D::Solve(), one or both of "
         "`x` and `y` were empty!");

    // initialize f[x][y] to zeros with proper dimensions
    std::vector< std::vector<T> > f( x.size(), std::vector<T>(y.size(), 0.0));

    // omp_set_num_threads() should be called prior to here!
    #pragma omp parallel for shared(f)
    for (std::size_t i = 0; i < x.size(); i++){

		if ( verbose > 2 && !omp_get_thread_num() )
            display -> Progress(i, x.size(), omp_get_num_threads() );

		for (std::size_t j = 0; j < y.size(); j++){

	        T sum = 0.0;

	        for (std::size_t k = 0; k < N; k++){

	            T WW     = W(x[i] - _x[k], y[j] - _y[k]);
	            f[i][j] += WW * _z[k];
	            sum     += WW;
	        }

			if (unbiased){

				// adjust the sum over weights to `unbias` the result
				// only relavent when called from Stdev()!!!
				f[i][j] /= (1.0 - 1.0 / N) * sum;

			} else f[i][j] /= sum;
	    }
	}

    return f;
}

template<class T>
std::vector< std::vector<T> > KernelFit2D<T>::Variance(const std::vector<T> &x,
	const std::vector<T> &y, const bool unbiased){

	//
	// Solve for the estimated standard deviation by evaluating
	// the profile *at* the raw data points.
	//

	if ( x.empty() || y.empty() )
		throw KernelFitError("From KernelFit2D::Variance(), one or both of the "
	    "input vectors were empty!");

	// initialize vector for profile at data points
	std::vector<T> f(N, 0.0);

	// solve profile at data points
	#pragma omp parallel for shared(f)
	for (std::size_t i = 0; i < N; i++){

		if ( verbose > 2 && !omp_get_thread_num() )
            display -> Progress(i, N, omp_get_num_threads() );

	    T sum = 0.0;

	    for (std::size_t j = 0; j < N; j++){

	        T W   = Kernel(_x[i] - _x[j], _y[i] - _y[j]);
	        f[i] += W * _z[j];
	        sum  += W;
	    }

	    f[i] /= sum;
	}

	if (verbose > 2)
		display -> Progress(1, 1); // complete

	// solve for variances at data points
	std::vector<T> var(N, 0.0);
	for (std::size_t i = 0; i < N; i++)
		var[i] = std::pow(_z[i] - f[i], 2.0);

	// solve for smooth surface through variance points
	KernelFit2D<T> profile(_x, _y, var, _b);

	return profile.Solve(x, y, unbiased);
}

template<class T>
std::vector< std::vector<T> > KernelFit2D<T>::Variance(const std::vector<T> &x,
	const std::vector<T> &y, T (*W)(T, T), const bool unbiased){

	//
	// Solve for the estimated standard deviation by evaluating
	// the profile *at* the raw data points. This version accepts an
	// alternative kernel function provided by the user.
	//

	if ( x.empty() || y.empty() )
		throw KernelFitError("From KernelFit2D::Variance(), one or both of the "
	    "input vectors were empty!");

	// initialize vector for profile at data points
	std::vector<T> f(N, 0.0);

	// solve profile at data points
	#pragma omp parallel for shared(f)
	for (std::size_t i = 0; i < N; i++){

		if ( verbose > 2 && !omp_get_thread_num() )
            display -> Progress(i, N, omp_get_num_threads() );

	    T sum = 0.0;

	    for (std::size_t j = 0; j < N; j++){

	        T WW  = W(_x[i] - _x[j], _y[i] - _y[j]);
	        f[i] += WW * _z[j];
	        sum  += WW;
	    }

	    f[i] /= sum;
	}

	if (verbose > 2)
		display -> Progress(1, 1); // complete
	
	// solve for variances at data points
	std::vector<T> var(_x.size(), 0.0);
	for (std::size_t i = 0; i < _x.size(); i++)
		var[i] = std::pow(_z[i] - f[i], 2.0);

	// solve for smooth surface through variance points
	KernelFit2D<T> profile(_x, _y, var, _b);

	return profile.Solve(x, y, unbiased);
}

template<class T>
std::vector< std::vector<T> > KernelFit2D<T>::StdDev(const std::vector<T> &x,
	const std::vector<T> &y, const bool unbiased){

	//
	// Solve for the estimated standard deviation by evaluating
	// the profile *at* the raw data points.
	//

	if ( x.empty() || y.empty() )
		throw KernelFitError("From KernelFit2D::StdDev(), one or both of the "
	    "input vectors were empty!");

	// initialize vector for profile at data points
	std::vector< std::vector<T> > stdev = Variance(x, y, unbiased);

	// take the sqrt for the standard deviation
	#pragma omp parallel for shared(stdev)
	for (std::size_t i = 0; i < x.size(); i++)
	for (std::size_t j = 0; j < y.size(); j++)
		stdev[i][j] = std::sqrt(stdev[i][j]);
}

template<class T>
std::vector< std::vector<T> > KernelFit2D<T>::StdDev(const std::vector<T> &x,
	const std::vector<T> &y, T (*W)(T, T), const bool unbiased){

	//
	// Solve for the estimated standard deviation by evaluating
	// the profile *at* the raw data points. This version accepts an
	// alternative kernel function provided by the user.
	//

	if ( x.empty() || y.empty() )
		throw KernelFitError("From KernelFit2D::StdDev(), one or both of the "
	    "input vectors were empty!");

	// initialize vector for profile at data points
	std::vector< std::vector<T> > stdev = Variance(x, y, W, unbiased);

	// take the sqrt for the standard deviation
	#pragma omp parallel for shared(stdev)
	for (std::size_t i = 0; i < x.size(); i++)
	for (std::size_t j = 0; j < y.size(); j++)
		stdev[i][j] = std::sqrt(stdev[i][j]);

	return stdev;
}

// template classes
template class KernelFit1D<float>;
template class KernelFit2D<float>;
template class KernelFit1D<double>;
template class KernelFit2D<double>;
template class KernelFit1D<long double>;
template class KernelFit2D<long double>;

} // namespace Gaia
