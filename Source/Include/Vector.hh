// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/Vector.hh
//
// This is the header file for the `Vector` class. A `Vector` in this case
// is the simple representation of a point in R3. This allows the `Profile`
// classes to have a uniform implimentation for taking position information.
// `Theta` here is the polar angle, and `Phi` is the polar angle.

#ifndef _VECTOR_HH_
#define _VECTOR_HH_

namespace GAIA {

class Vector {

public:

	Vector(){}
	Vector(double x = 0.0, double y = 0.0, double z = 0.0);

	double X()     const; // return x coord.
	double Y()     const; // return y coord.
	double Z()     const; // return z coord.
	double R()     const; // return radius (cylindrical!)
	double Rho()   const; // return radius (spherical!)
	double Phi()   const; // return azimuthal angle (radians)
	double Theta() const; // return polar angle (radians)
	
	// set coords.
	void SetX(double x);
	void SetY(double y);
	void SetZ(double z);
	void SetXYZ(double x, double y, double z);
	
protected:
	
	// coordinates
	double _x, _y, _z;

};

} // namespace GAIA

#endif