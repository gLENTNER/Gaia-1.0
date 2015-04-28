// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/Profiles.hh
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

#include <ProfileBase.hh>
#include <Vector.hh>

namespace Gaia {

// CREATE USER DEFINED PROFILES HERE (INSIDE THE `Gaia` NAMESPACE)  ...
// -----------------------------------------------------------------------

// Profile for the mass density of the galaxy disk (isotropic)
class MassDensity: public ProfileBase {
public:

    MassDensity(): ProfileBase("MassDensity"){ }

    virtual double Function(const Vector& p){

        double no     = 0.500; // normalization coefficient
        double Rs     = 4.250; // scale radius for the disk
        double Zthick = 0.750; // scale height of the thick disk
        double Zthin  = 0.015; // scale height of the thin disk
        double m      = 0.020; // relative scale, thick disk

        return no * exp( -pow(p.R(), 2.0) / (2.0 * Rs * Rs) ) * (
            exp( -pow(p.Z(), 2.0) / (2.0 * Zthin * Zthin) ) +
            m * exp( -pow(p.Z(), 2.0) / (2.0 * Zthick * Zthick) ) );
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
