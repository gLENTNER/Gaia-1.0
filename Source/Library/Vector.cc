// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/Vector.hh
//
// This is the source file for the `Vector` class. A `Vector` in this case
// is the simple representation of a point in R3. This allows the `Profile`
// classes to have a uniform implimentation for taking position information.
// `Theta` here is the polar angle, and `Phi` is the polar angle.

#include <cmath>
#include <Vector.hh>
#include <Exception.hh>

#ifndef Half_Pi
// no need for repeat calculations
#define Half_Pi 1.5707963267948966
#endif

namespace GAIA {

// simple accessors, `getters`
inline double Vector::X() const {return _x;}
inline double Vector::Y() const {return _y;}
inline double Vector::Z() const {return _z;}

// return the radius (cylindrical)
inline double Vector::R() const {return sqrt(_x*_x + _y*_y);}

// return the radius (spherical) and/or magnitude
inline double Vector::Rho() const {return sqrt(_x*_x + _y*_y + _z*_z);}

// return the azimuthal angle
inline double Vector::Phi() const {
	if      (x != 0) return atan2(_y, _x);
	// in this case, it's easy to work out what the angle should be
	else if (y == 0) return 0.0;
	else if (y  > 0) return  Half_Pi;
	else             return -Half_Pi;
}

// return the polar angle
inline double Vector::Theta() const {
	if (_x != 0.0 || _y != 0 || _z != 0) return acos(_z / R());
	// in this program, R() should never be exactly zero;
	// the odds of getting identically 0.0 for all x, y, and z is
	// astronomically low. 
	// FIXME: find an alternative for Theta()
	else throw DivError("Radius was zero in Vector::Theta()!");
}

// `setters` (cartesian coordinates only)
inline void Vector::SetX(const double& x){_x = x;}
inline void Vector::SetY(const double& y){_y = y;}
inline void Vector::SetZ(const double& z){_z = z;}
inline void Vector::SetXYZ(const double& x, const double&y, const double& z){
	_x = x; 
	_y = y; 
	_z = z;
}

} // namespace GAIA