// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Library/ProfileBase.cc
//
// #TODO:30 source

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include <ProfileBase.hpp>
#include <Exception.hpp>
#include <Interpolate.hpp>
#include <Parser.hpp>

namespace Gaia {

// construct the function maps
ProfileBase::ProfileBase(const std::string &name,
	const std::string &axis1, const std::string &axis2){

	// set pointers to nullptr immediately
	Linear_Data   = nullptr;
	BiLinear_Data = nullptr;

	// set name
	_name = name;

	// build coordinate function map
	Coord["X"] = X;
	Coord["Y"] = Y;
	Coord["Z"] = Z;
	Coord["R"] = R;
	Coord["Rho"] = Rho;
	Coord["Phi"] = Phi;
	Coord["Theta"] = Theta;

	// specified axis
	_axis1 = axis1;
	_axis2 = axis2;

	// grab the parser
	Parser *parser = Parser::GetInstance();

	// map for choosing limits dynamically
	Limits["X"] = parser -> GetXlimits();
	Limits["Y"] = parser -> GetYlimits();
	Limits["Z"] = parser -> GetZlimits();

	// assume analytical
	_analytical = true;
}

ProfileBase::~ProfileBase(){

	if (Linear_Data){
		delete Linear_Data;
		Linear_Data = nullptr;
	}

	if (BiLinear_Data){
		delete BiLinear_Data;
		BiLinear_Data = nullptr;
	}
}

void ProfileBase::Initialize(std::string &filename){

	//
	// read in contents of file and parse data
	//

	// we are not going to use the Function() so not analytical
	_analytical = false;

	// grab variables from parser
	Parser *parser = Parser::GetInstance();
	int verbose = parser -> GetVerbosity();

	// update the user
	if (verbose) std::cout
		<< " Initializing `" << _name << "` profile from file `"
		<< filename << "` ...";

    ReplaceAll("~", std::string(getenv("HOME")), filename);
	std::ifstream input( filename.c_str() );

	if ( input ) {

		for ( std::string line; std::getline(input, line); )
			_data.push_back( ReadElements(line) );

	} else throw IOError("`"+filename+"` failed to open properly!\n");

	// ensure appropriate input (dimensionally)
	for ( std::size_t i = 1; i < _data.size(); i++ ){

		// all rows must have the same length!
		if ( _data[i].size() != _data[i-1].size() ){

			std::stringstream warning;
			warning << "From file `" << filename << "`, rows " << i - 1;
			warning << " and " << i << " don't have the same number of ";
			warning << " elements!\n";

			throw ProfileError( warning.str() );
		}
	}

	if ( _data.size() < 2 || _data[0].size() < 2 ){

		std::stringstream warning;
		warning << "From file `" << filename << "`, there must be ";
		warning << "at least two rows and two columns. See README.md ";
		warning << "file for more information on how to create your ";
		warning << "input files!\n";

		throw ProfileError( warning.str() );
	}

	if ( _data.size() == 2 || _data[0].size() == 2 ){

		// we are a 1D profile
		_1D = true;
		_2D = false;

		if (_data.size() == 2){

			// horizontal
			_x = _data[0];
			_y = _data[1];

		} else {

			// vertical
			for ( std::size_t i = 0; i < _data.size(); i++ ){
				_x.push_back( _data[i][0] );
				_y.push_back( _data[i][1] );
			}
		}

		// check constructor arguments
		if ( _axis1.empty() ){

			std::stringstream warning;
			warning << "From file `" << filename << "`, I have detected ";
			warning << " a 1D data set but you didn't specify an axis ";
			warning << " for `" << _name << "` in `Profiles.hpp`!\n";

			throw ProfileError( warning.str() );

		} else if ( Coord.find(_axis1) == Coord.end() ){

			std::stringstream warning;
			warning << "In `" << _name << "` from `Profiles.hpp`, ";
			warning << "the axis specified does not match any of the ";
			warning << "available coordinates!\n";

			throw ProfileError( warning.str() );

		} else if ( !(_axis2.empty()) ){

			std::stringstream warning;
			warning << "From file `" << filename << "`, I have detected ";
			warning << " a 1D data set but you specified two axis ";
			warning << " for `" << _name << "` in `Profiles.hpp`! ";
			warning << "This is ambiguous, see README.md for details.\n";

			throw ProfileError( warning.str() );
		}

		// construct linear data member
		Linear_Data = new Interpolate::Linear<double>(_x, _y);

	} else {

		// we are a two dimensional surface profile
		_1D = false;
		_2D = true;

		// check constructor axis arguments
		if ( _axis1.empty() || _axis2.empty() ){

			std::stringstream warning;
			warning << "From file `" << filename << "`, I have detected ";
			warning << "a 2D data set but you didn't specify two axis ";
			warning << "for `" << _name << "` in `Profiles.hpp`!\n";

			throw ProfileError( warning.str() );

		} else if ( Coord.find(_axis1) == Coord.end() ){

			std::stringstream warning;
			warning << "In `" << _name << "` from `Profiles.hpp`, ";
			warning << "the first axis specified does not match any ";
			warning << "of the available coordinates!\n";

			throw ProfileError( warning.str() );

		} else if ( Coord.find(_axis2) == Coord.end() ){

			std::stringstream warning;
			warning << "In `" << _name << "` from `Profiles.hpp`, ";
			warning << "the second axis specified does not match any ";
			warning << "of the available coordinates!\n";

			throw ProfileError( warning.str() );

        } else if ( Limits.find(_axis1) == Limits.end() ||
                   Limits.find(_axis2) == Limits.end() ){

            std::stringstream warning;
            warning << "In `" << _name << "` from `Profiles.hpp`, one or ";
            warning << "more of the axis specified does not match any ";
            warning << "of the allowed coordinates!\n";

            throw ProfileError( warning.str() );
        }

        // the `x` and `y` are now line-spaces
        _x = Linespace(Limits[_axis1][0], Limits[_axis1][1], _data[0].size());
        _y = Linespace(Limits[_axis2][0], Limits[_axis2][1], _data.size());

        // construct the 2D interpolation object
        BiLinear_Data = new Interpolate::BiLinear<double>(_x, _y, _data);
	}

	// update user
	if (verbose){
		std::string dim = _1D ? " (1D) " : " (2D) ";
		std::cout << dim << _data.size() << " x " << _data[0].size() << "\n";
	}
}

// generalized function for evaluating the profile
double ProfileBase::Evaluate(const Vector &vec){

		if (_analytical) return Function(vec);

		if (_1D) return Linear_Data -> Interpolate( Coord[_axis1](vec) );

		else return BiLinear_Data -> Interpolate( Coord[_axis1](vec),
						Coord[_axis2](vec) );
}

// convert new line of text from file into vector<double>
std::vector<double> ProfileBase::ReadElements(std::string &line){

	// replace all `,` delimiters (if any) to spaces
	ReplaceAll(",", " ", line);

	// stream reads elements from string
	std::stringstream buffer(line);

	// empty vector to be filled
	std::vector<double> data;

	for ( double element; buffer >> element; )
		data.push_back(element);

	return data;
}

std::vector<double> ProfileBase::Linespace(const double start, const double end,
	const std::size_t length){

	// initialize the vector
	std::vector<double> linespace( length, 0.0 );

	linespace[0] = start;
	double dx    = (end - start) / double(length - 1);

	for ( std::size_t i = 1; i < length; i++ )
		linespace[i] = linespace[i-1] + dx;

	return linespace;
}

// replace all instances of `search_str` in `input_str` with `replace_str`
void ProfileBase::ReplaceAll(const std::string &search_str,
    const std::string &replace_str, std::string& input_str ){

    std::size_t pos = 0;

    while ( ( pos = input_str.find(search_str, pos) ) != std::string::npos ){

        input_str.replace( pos, search_str.length( ), replace_str );
        pos += replace_str.length( );
    }
}

} // namespace Gaia
