// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/ProfileBase.hpp
//
// TODO: header

#ifndef _PROFILEBASE_HH_
#define _PROFILEBASE_HH_

#include <string>
#include <vector>
#include <map>

#include <Vector.hpp>
#include <Interpolate.hpp>
#include <Parser.hpp>

namespace Gaia {

class ProfileBase {

public:

	// construct `Profile` with a name and axis information
	ProfileBase(const std::string &name, const std::string &axis1 = "",
		const std::string &axis2 = "");
	
    ~ProfileBase();

	void Initialize(std::string &filename);

	std::string Name(){return _name;}

	// analytical `Function` is unique to each derived `Profile`
	virtual double Function(const Vector &vec){ return 1.0; }

	// generalized accessor function chooses what to do
	double Evaluate(const Vector &vec);

private:

	// convert string to vector<double>
	std::vector<double> ReadElements(std::string &line);
    // replace elements of one string in another
    void ReplaceAll(const std::string&, const std::string&, std::string&);

	// used to differentiate derived `Profile` classes at runtime
	std::string _name;

	// data from file
	std::vector< std::vector<double> > _data;

	// 1D data (`x` is not necessarily x and y = f(x) )
	std::vector<double> _x, _y;

	// limits needed to construct line-space given a 2D profile
	std::map< std::string, std::vector<double> > Limits;

	// map of functions, axis1 and axis2
	std::map< std::string, double (*)(const Vector&) > Coord;

	// functions in the above map (calls to vector coordinates)
	static double X(const Vector &vec)     { return vec.X();     }
	static double Y(const Vector &vec)     { return vec.Y();     }
	static double Z(const Vector &vec)     { return vec.Z();     }
	static double R(const Vector &vec)     { return vec.R();     }
	static double Rho(const Vector &vec)   { return vec.Rho();   }
	static double Phi(const Vector &vec)   { return vec.Phi();   }
	static double Theta(const Vector &vec) { return vec.Theta(); }

	// given axis from derived class constructor
	std::string _axis1, _axis2;

	// flags to signify 1D or 2D
	bool _1D, _2D;

	// flag to signify analyticity
	bool _analytical;

	// member interpolators
	Interpolate::Linear<double>   *Linear_Data;
	Interpolate::BiLinear<double> *BiLinear_Data;

	// helper function for Initialize()
	std::vector<double> Linespace(const double, const double, const std::size_t);
};

} // namespace Gaia

#endif
