// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/ProfileBase.hh
//
// TODO: header

#ifndef _PROFILEBASE_HH_
#define _PROFILEBASE_HH_

#include <string>
#include <vector>

#include <Vector.hh>

namespace GAIA {

class ProfileBase {

public:
	
	// construct `Profile` with a name
	ProfileBase(const std::string &name);
	~ProfileBase();
	
	void Initialize(const std::string &filename);
	
	void Evaluate(const Vector&);
	
	std::string Name(){return _name;}
	
	// analytical `Function` is unique to each derived `Profile`
	virtual double Function(const Vector &position){return 1;}
	
private:
	
	// convert string to vector<double>
	std::vector<double> ReadElements(std::string &line);
	
	// used to differentiate derived `Profile` classes at runtime
	std::string _name;
	
	// data from file
	std::vector< std::vector<double> > _data;
	
	// 1D data (`x` is not necessarily x and y = f(x) )
	std::vector<double> _x, _y;
	
	// map of functions, dimensions
	std::map< std::string, double (*)(const Vector&) > _coords;
	
	// functions in the above map (calls to vector coordinates)
	double X(const Vector &coord)     { return coord.X();     }
	double Y(const Vector &coord)     { return coord.Y();     }
	double Z(const Vector &coord)     { return coord.Z();     }
	double R(const Vector &coord)     { return coord.R();     }
	double Rho(const Vector &coord)   { return coord.Rho();   }
	double Phi(const Vector &coord)   { return coord.Phi();   }
	double Theta(const Vector &coord) { return coord.Theta(); }
	
};

} // namespace GAIA

#endif
