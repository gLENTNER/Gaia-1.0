// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Library/ProfileBase.cc
//
// TODO: source

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <ProfileBase.hh>
#include <Strings.hh>
#include <Exception.hh>

namespace GAIA {

// construct the function maps
ProfileBase::ProfileBase(const std::string &name, 
	const std::string &dim1 = "", const std::string &dim2 = ""){
	
	// set name
	_name = name;
	
	// build coordinate function map
	_dimension["X"] = X;
	_dimension["Y"] = Y;
	_dimension["Z"] = Z;
	_dimension["R"] = R;
	_dimension["Rho"] = Rho;
	_dimension["Phi"] = Phi;
	_dimension["Theta"] = Theta;
	
	// grab the parser
	Parser *parser = Parser::GetInstance();
	
	// map for choosing limits dynamically
	_limits["X"] = parser -> GetXlimits;
	_limits["Y"] = parser -> GetYlimits;
	_limits["Z"] = parser -> GetZlimits;
}

void ProfileBase::Initialize(const std::string &filename,
	const std::string &dim1 = "", const std::string &dim2 = ""){
	
	// read in contents of file and send to appropriate function
	// for parsing ...
	
	// grab variables from parser
	Parser *parser = Parser::GetInstance();
	int verbose = parser -> GetVerbosity();
	
	// update the user
	if (verbose) std::cout
		<< " Initializing `" << _name << "` profile from file `"
		<< filename << "` ...";
	
	std::ifstream input( filename.c_str() );
	
	if ( input ) {
		
		for ( std::string line; std::getline(input, line); )
			_data.push_back( ReadElements(line) );		
		
	} else throw IOError("`" + filename + "` failed to open properly!");
	
	// ensure appropriate input (dimensionally)
	for ( int i = 1; i < _data.size(); i++ ){
		
		// all rows must have the same length!
		if ( _data[i].size() != _data[i-1].size() ){
			
			std::stringstream warning;
			warning << "From file `" << filename << "`, rows " << i - 1;
			warning << " and " << i << " don't have the same number of ";
			warning << " elements!";
			
			throw InputError( warning.str() );
		}
	}
	
	if ( _data.size() == 2 || _data[0].size() == 2 ){
		// we are a 1D profile
		_oneD = true; 
		_twoD = false;
		if (_data.size() == 2){
			// horizontal
			_x = _data[0];
			_y = _data[1];
		
		} else {
			// vertical
			for ( int i = 0; i < _data.size(); i++){
				_x.push_back( _data[i][0] );
				_y.push_back( _data[i][1] );
			}
		}
	
	} else {
		// we are a two dimensional surface profile
		_oneD = false;
		_twoD = true;
		// the `x` and `y` are now line-spaces
		_dimension[1]
		std::vector<double> xlimits = _limits[dim1];
		std::vector<double> ylimits = _limits[dim2];
		_x = Linespace(xlimits[0], xlimits[1], _data.size());
		_y = Linespace(ylimits[0], ylimits[1], _data[0].size());
	}
	
	// update user
	if (verbose){
		std::string dim = _oneD ? " (1D) " : " (2D) ";
		std::cout << dim << _data.size() << " x " << _data[0].size() << "\n";
	}
}

// convert new line of text from file into vector<double>
std::vector<double> ProfileBase::ReadElements(std::string &line){
	
	// replace all `,` delimiters (if any) to spaces
	ReplaceAll(",", " ", line);
	
	// stream reads elements from string
	std::stringstream input(line);
	
	// empty vector to be filled
	std::vector<double> data;
	
	for ( double element; input >> element; )
		data.push_back(element);
			
	return data;
}

} // namespace GAIA
