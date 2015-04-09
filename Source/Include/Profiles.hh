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

namespace GAIA {
	
// CREATE USER DEFINED PROFILES HERE (INSIDE THE `GAIA` NAMESPACE)  ...
// -----------------------------------------------------------------------

// Profile for the mass density of the galaxy disk (isotropic in the disk)
class MassDensity: public ProfileBase {
public:
	
	MassDensity(): ProfileBase("MassDensity"){ }
	
	virtual double Function(const Vector& p){
	
		double N_0    = 0.01; // normalization coefficient
		double R_D    = 6000; // scale radius for the disk
		double Zthick = 1000; // scale height of the thick disk
		double Zthin  = 300;  // scale height of the thin disk
	
		return N_0 * exp( -p.R() / R_D ) * ( exp( -p.Z() / Zthick ) + 
			0.02 * exp( -p.Z() / Zthin ) );
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

// Empirical profile read in from file 
class Metallicity : public ProfileBase {
public:
	
	Metallicity(): ProfileBase("Metallicity", "R"){ }
};

// Example for NGC1300 from HST FITS image ...
class SurfaceProfile: public ProfileBase {
public:
	
	SurfaceProfile(): ProfileBase("SurfaceProfile", "X", "Y"){ }
};

} // namespace GAIA