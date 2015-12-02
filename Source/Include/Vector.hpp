// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// GNU General Public License v3.0
// Include/Vector.hpp
//
// This is the header file for the `Vector` class. A `Vector` in this case
// is the simple representation of a point in R3. This allows the `Profile`
// classes to have a uniform implimentation for taking position information.
// `Theta` here is the polar angle, and `Phi` is the polar angle. The entire
// class is defined in the header for efficiency.

#ifndef _VECTOR_HH_
#define _VECTOR_HH_

#include <cmath>
#include <iostream>

#include <Exception.hpp>

#ifndef Half_Pi
// no need for repeat calculations
#define Half_Pi 1.5707963267948966
#endif

namespace Gaia {

class Vector {

public:

	Vector(double x = 0.0, double y = 0.0, double z = 0.0):
		_x(x), _y(y), _z(z){ }

	// simple accessors, `getters`
	double X() const {return _x;}
	double Y() const {return _y;}
	double Z() const {return _z;}

	// return the radius (cylindrical)
	double R() const {return sqrt(_x*_x + _y*_y);}

	// return the radius (spherical) and/or magnitude
	double Rho() const {return sqrt(_x*_x + _y*_y + _z*_z);}

	// return the azimuthal angle
	double Phi() const {
		if (_x != 0) return atan2(_y, _x);
		// in this case, it's easy to work out what the angle should be
		else if (_y == 0) return 0.0;
		else if (_y  > 0) return  Half_Pi;
		else              return -Half_Pi;
	}

	// return the polar angle
	double Theta() const {
		if (_x != 0.0 || _y != 0.0 || _z != 0.0) return acos(_z / Rho());
		// in this program, R() should never be exactly zero;
		// the odds of getting identically 0.0 for all x, y, and z is
		// astronomically low.
		// FIXME: find an alternative for Theta()
		else throw DivError("Radius was zero in Vector::Theta()!");
	}

	// magnitude is simply an alias for `Rho`
	double Mag() const {return sqrt(_x*_x + _y*_y + _z*_z);}

	// `setters` (cartesian coordinates only)
	void SetX(const double& x){_x = x;}
	void SetY(const double& y){_y = y;}
	void SetZ(const double& z){_z = z;}
	void SetXYZ(const double& x, const double&y, const double& z){
		_x = x; _y = y; _z = z;
	}

	// addition
	Vector operator+ ( const Vector &rhs ){

	        Vector result( X() + rhs.X(), Y() + rhs.Y(), Z() + rhs.Z() );
	        return result;
	    }

	// subtraction
	Vector operator- ( const Vector &rhs ){

		Vector result( X() - rhs.X(), Y() - rhs.Y(), Z() - rhs.Z() );
		return result;
	}

	// output to a stream
	friend std::ostream& operator<< (std::ostream &aStream, const Vector &aVec){
		aStream << aVec.X() << " " << aVec.Y() << " " << aVec.Z();
		return aStream;
	}

protected:

	// coordinates
	double _x, _y, _z;

};

} // namespace Gaia

#endif
