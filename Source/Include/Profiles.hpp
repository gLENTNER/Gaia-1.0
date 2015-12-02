// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// GNU General Public License v3.0
// Include/Profiles.hpp
//
// The various `Profiles` can be defined here. Every `Profile` class should
// be derived from `ProfileBase`. If your profile will import from a file there
// is nothing to do but declare it. If your profile has some unique analytical
// form, overload the `Function` definition. The template is
//
//     virtual double Function(const Vector&);
//
// Following and/or remake the below examples... (and stay in the namespace!).

#include <cmath>

#include <ProfileBase.hpp>
#include <Vector.hpp>

namespace Gaia {

// CREATE USER DEFINED PROFILES HERE (INSIDE THE `Gaia` NAMESPACE)  ...
// -----------------------------------------------------------------------

// Profile for the mass density of the galactic disk (isotropic),
// Parameterization taken from McMillan, 2011.
class MilkyWay: public ProfileBase {
public:

    MilkyWay(): ProfileBase("MilkyWay"){ }

    virtual double Function(const Vector& p){

        // constant in front is pseudo-normalization parameter
        return 0.2 * (

            // thin disk parameterization:
            // z_d = 0.3 +/- 0.?? kpc
            // r_d = 2.6 +/- 0.52 kpc
            std::exp( -std::abs( p.Z() ) / 0.3 - p.R() / 2.6 ) / 0.3

            +

            // thick disk parameterization:
            // z_d = 0.9 +/- 0.?? kpc
            // r_d = 3.6 +/- 0.72 kpc
            std::exp( -std::abs( p.Z() ) / 0.9 - p.R() / 3.6 ) / 0.9
            );
    }
};

// model spirals
class Spiral: public ProfileBase {
public:

    Spiral(): ProfileBase("Spiral"){ }

    virtual double Function(const Vector &p){

        double n   = 1.0;
        double Rs  = 16.863;
        double xi  =  1.5;
        double pi2 =  3.141592653589793 * 2.0;
        double tau =  2.0;

        return pow(cos(n * p.Phi() - xi * pi2 * p.R() / Rs), tau);
    }
};

class Metallicity: public ProfileBase {
public:

    Metallicity(): ProfileBase("Metallicity"){ }

    virtual double Function( const Vector& p ){

        double No    =  0.452322; // normalization coefficient
        double base  =  0.760000; // base level
        double slope =  0.880000; // over all slope
        double Rgal  = 18.759116; // semi-major axis of galaxy
        double Rmid  =  9.379558; // center of peak metallicity
        double Rs    =  3.126519; // scale radius

        return No * (base + slope * p.R() / Rgal +
            exp( -pow(p.R() - Rmid, 2.0) / (2 * Rs * Rs)));
    }
};

// Hypothetical profile for habitability in the disk (radial)
class Habitability: public ProfileBase {
public:

	Habitability(): ProfileBase("Habitability"){ }

	virtual double Function(const Vector& position){

		double N_0   = 0.01; // normalization coefficient
		double sigma = 300;  // bandwidth for profile
		double R_c   = 7500; // orbit of co-rotation

		// Gaussian around orbit of co-rotation
		return N_0 * exp(-pow(position.R() - R_c, 2.0) / (2.*sigma*sigma));
	}
};

// Example for NGC1300 from HST FITS image ...
class Surface: public ProfileBase {
public:

	Surface(): ProfileBase("Surface", "X", "Y"){ }
};

} // namespace Gaia
