// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// GNU General Public License v3.0
// Include/Interpolation.hpp
//
// General interpolation objects. Interp1D and Interp2D provide
// Linear and Bilinear Interpolation respectively.
//
// These object throw an InterpError exception derived from the
// std::exception

#ifndef _INTERPOLATE_HH_
#define _INTERPOLATE_HH_

#include <exception>
#include <string>
#include <vector>

namespace Gaia {

namespace Interpolate {

template<class T>
class Linear {

public:

    Linear(const std::vector<T>& x, const std::vector<T> &y);

    // find a new `y` vector given a new `x` vector
    std::vector<T> Interpolate(const std::vector<T> &x);

    // find a new `y` given a new `x`
    T Interpolate(const T &x);

private:

    // keep local data
    std::vector<T> x, y, m;

};

template<class T>
class BiLinear {

public:

    BiLinear(const std::vector<T>& x, const std::vector<T> &y,
        const std::vector< std::vector<T> > &z);

    // find a new `z` matrix given a new `x` and `y` grid
    std::vector< std::vector<T> > Interpolate(const std::vector<T> &x,
        const std::vector<T> &y);

    // find a new `z` given a new `x`, `y` pair
    T Interpolate(const T &x, const T &y);

private:

    // set of two vectors represent rectilinear grid
    std::vector<T> x, y;

    // matrix of values for x, y pairs
    std::vector< std::vector<T> > z;

    // set up 1D linear interpolation objects for each row
    std::vector< Linear<T> > partial;

};

// base exception class for Interpolator objects
class InterpException : public std::exception {

public:

	explicit InterpException(const std::string& msg): _msg(msg){ }
	virtual ~InterpException() throw() { }
	virtual const char* what() const throw(){ return _msg.c_str(); }

protected:

	std::string _msg;
};

// exception thrown by Interpolation objects
class InterpError : public InterpException {
public:

    InterpError(const std::string& msg): InterpException(
		"\n --> InterpError: " + msg){ }
};

} // namespace Interpolate

} // namespace Gaia

#endif
