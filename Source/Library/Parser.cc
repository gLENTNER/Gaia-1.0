// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Library/Parser.cc
//
// This source file contains the method definitions for the `Parser` object.
// The purpose of this singleton class is to interpret the arguments passed
// to `main` at runtime and to read the configuration file, retaining the
// parameters for retrieval by the other objects.
//
// TODO: Add `source` functionality for rc files, _random_seed

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <list>
#include <omp.h>
#include <stdlib.h>

#include <Parser.hh>
#include <Exception.hh>

namespace Gaia {

// static pointer for singleton class
Parser* Parser::instance = NULL;

// retrieval of pointer
Parser* Parser::GetInstance() {

	if ( !instance )
		instance = new Parser();

	return instance;
}

// destroy the singleton
void Parser::Release( ) {

	if ( instance ) {
		delete instance;
		instance = NULL;
	}
}

// set up the simulation environment
void Parser::Setup(const int argc, const char *argv[]){

	// display usage
	if ( argc == 1 ) throw Usage(
    "Gaia [--num-particles=] [--num-trials=] [--num-threads=] [--set-verbose=0|1|2|3]\n"
    "\t[--out-path=] [--raw-path=] [--tmp-path=] [--pos-path=] [--first-seed=]\n"
    "\t[--rc-file=] [--no-analysis] [--keep-raw] [--keep-pos] [--debug] \n\n\t"
    "An application for building 3D numerical models of systems of particles\n\t"
    "using a Monte Carlo rejection chain algorithm based on probability density\n\t"
    "functions (PDFs) defined by the user. A nearest neighbor analysis is \n\t"
    "performed on each of a number of trial constructions for the system.\n\n\t"
    "See the README.md file for more detailed usage information.\n");

	// turn argv into vector
	for (int i = 1; i < argc; i++)
		_cmd_args.push_back( std::string(argv[i]) );

	// define `Command` map
	// Command["set"]     = Set;
	// Command["include"] = Include;

	// set the defaults for the parameters
	SetDefaults();

	// read configuration file
	ReadRC();

	// allow reassignment from command line (done after rc file)
	Interpret();

	// check arguments
	Rectify();
}


// set the defaults for the parameters
void Parser::SetDefaults(){

	// default arguments
	argument["--num-particles"] = "~";   // necessarily reset
	argument["--num-trials"   ] = "30";
	argument["--num-threads"  ] = "1";
	argument["--set-verbose"  ] = "2";
	argument["--out-path"     ] = "Gaia-out-";
	argument["--raw-path"     ] = "Gaia-raw-";
	argument["--tmp-path"     ] = "Gaia-tmp-";
	argument["--pos-path"     ] = "Gaia-pos-";
	argument["--no-analysis"  ] = "0";
	argument["--keep-raw"     ] = "0";
	argument["--keep-pos"     ] = "0";
	argument["--rc-file"      ] = "~/.Gaiarc";
	argument["--first-seed"   ] = "~"; // will be assigned if not given
	argument["--debug"        ] = "0";

	// arguments who don't need an assigment
	implicit["--no-analysis"] = "~";
	implicit["--keep-raw"   ] = "~";
	implicit["--keep-pos"   ] = "~";
	implicit["--debug"      ] = "~";

	// list values as `not given` before assignments
	_given_xlims = _given_ylims  = _given_zlims = false;
	_given_xy    = _given_radial = false;
}

// parse the commands in the configuration (rc) file
void Parser::ReadRC(){

	// default location for `rc file`
	_rc_file = argument["--rc-file"];

	// look for alternative path to rc file in commandline arguments
	for ( std::vector<std::string>::iterator iter = _cmd_args.begin();
		iter != _cmd_args.end(); iter++ ){

		std::string arg( *iter );

		if ( arg.find("--rc-file=") != std::string::npos )
			_rc_file = arg.substr(arg.find("=") + 1, arg.length());
	}

	// replace `~` with `$HOME`
	ReplaceAll("~", std::string(getenv("HOME")), _rc_file);

	// list of commands by line
	std::list< std::vector<std::string> > command_list;

	// open rc file
	std::ifstream rc_file( _rc_file.c_str() );

	if (rc_file){

		// retrieve contents of file

		for ( std::string new_line; std::getline(rc_file, new_line); ){

			// strip comments
			Clip(new_line, "#");

			// create vector and append to command list
			command_list.push_back( Split(new_line) );
		}

	} else throw IOError("Failed to open `"+_rc_file+"`.");

	// keep count of lines
	_line_number = 0;

	// parse the commands
	for (std::list< std::vector<std::string> >::iterator
		cmd = command_list.begin(); cmd != command_list.end(); cmd++){

		_line_number++;
		std::vector<std::string> line = *cmd;

		// skip empty lines
		if (line.size() == 0) continue;

		// check that the first word is appropriate
		//if ( Command.find(line[0]) == Command.end() ){
		if ( line[0] != "set" && line[0] != "include" ){

			std::stringstream warning;
			warning << "In file `" << _rc_file << "` on line " << _line_number;
			warning << ", `" << line[0] << "` was not a recognized ";
			warning << "command option!";
			throw InputError( warning.str() );
		}

		// give the command
		//Command[ line[0] ]( line, _rc_file, _line_number);

		if ( line[0] == "set" )
			Set(line);
		else
			Include(line);
	}
}

// parse the arguments passed to `main`
void Parser::Interpret() {

	// set `given` map to all false for all `arguments` not given in rc file
	for( std::map<std::string, std::string>::iterator arg = argument.begin();
		arg != argument.end(); arg++ )
		if( given.find(arg -> first) == given.end() )
			given[ arg -> first ] = false;

	// parse input arguments
	for ( std::vector<std::string>::iterator iter = _cmd_args.begin();
		iter != _cmd_args.end(); iter++ ){

		// string for this argument
		std::string arg( *iter );

		// split by '=' sign
		std::size_t pos = arg.find("=");
		if ( pos == std::string::npos ){
			// check if argument was an `implicit` argument
			if ( implicit.find(arg) != implicit.end() ){
				// the value is implicitely assumed
				given[arg] = true;

			} else {
				// this argument isn't understood
				throw InputError( "Missing assignment for " + arg + "!");
			}
		}

		// substring containing keyword argument
		std::string keyword = arg.substr(0, pos);

		// check for keyword argument in `argument` map
		if ( argument.find(keyword) == argument.end() )
			throw InputError( keyword + " is not a recognized parameter!" );

		// check that we have more slack in string
		if ( pos == arg.length( ) - 1 )
			throw InputError( "No assignment given for " + keyword + "!" );

		// value given by remainder of string
		std::string value = arg.substr( pos + 1, arg.length( ) - 1 );

		// valid assignment
		argument[keyword] = value;
		given[keyword]    = true;
	}
}

// rectify simulatin parameters
void Parser::Rectify(){

	// set `string` arguments
	_out_path = argument["--out-path"];
	_raw_path = argument["--raw-path"];
	_tmp_path = argument["--tmp-path"];
	_pos_path = argument["--pos-path"];
	_rc_file  = argument["--rc-file" ];

	// replace `~` with `$HOME`
	ReplaceAll("~", std::string(getenv("HOME")), _out_path);
	ReplaceAll("~", std::string(getenv("HOME")), _raw_path);
	ReplaceAll("~", std::string(getenv("HOME")), _tmp_path);
	ReplaceAll("~", std::string(getenv("HOME")), _pos_path);
	ReplaceAll("~", std::string(getenv("HOME")), _rc_file );

	// giving `--raw-path` or `--pos-path` implicitely means `--keep-*`
	_keep_raw = given["--raw-path"] || given["--keep-raw"] ? true : false;
	_keep_pos = given["--pos-path"] || given["--keep-pos"] ? true : false;
	_no_analysis = given["--no-analysis"] ? true : false;

	// `--no-analysis` makes no sense with these other options ...
	if ( _no_analysis && given["--out-path"] )
		throw InputError("With the `--no-analysis` flag, there will be no "
			"analysis conducted, but you specified an `--out-path`! ");
	if ( _no_analysis && _keep_raw )
		throw InputError("With the `--no-analysis` flag, there will be no "
			"analysis conducted, but you asked to keep `raw` files!");

	// stringstream used to convert between types
	std::stringstream convert;

	// check particle numbers
	if ( !given["--num-particles"] )
		throw InputError("User must provide --num-particles for system!");
	convert.str( argument["--num-particles"] );
	double as_double; // allows for scientific notation!
	if ( !( convert >> as_double ) || as_double < 2.0 )
		throw InputError("--num-particles must take integer value > 2!");
	_num_particles = as_double; // convert back to size_t

	// check verbosity
	convert.clear();
	convert.str( argument["--set-verbose"] );
	if ( !( convert >> _verbose ) || _verbose < 0 || _verbose > 3 )
		throw InputError("verbose takes 0, 1, 2, or 3.");

	// check num threads
	convert.clear();
	convert.str( argument["--num-threads"] );
	if ( !( convert >> _num_threads ) || _num_threads < 1 )
		throw InputError("--num-threads must be a positive integer!");
	if ( _num_threads > omp_get_max_threads() )
		throw InputError( "OpenMP says you have less than " +
			argument["--num-threads"] + " available!" );

	// check trial numbers
	convert.clear();
	convert.str( argument["--num-trials"] );
	if ( !( convert >> _num_trials ) || _num_trials < 1 )
		throw InputError("--num-trials needs a postive integer value!");

	// ensure we have Xlimits from rc file
	if ( !_given_xlims ) {
		std::stringstream warning;
		warning << "In file `" << _rc_file << "`, `Xlimits` was not given!";
		throw InputError( warning.str() );
	}

	// ensure we have Ylimits from rc file
	if ( !_given_ylims ) {
		std::stringstream warning;
		warning << "In file `" << _rc_file << "`, `Ylimits` was not given!";
		throw InputError( warning.str() );
	}

	// ensure we have Zlimits from rc file
	if ( !_given_zlims ) {
		std::stringstream warning;
		warning << "In file `" << _rc_file << "`, `Zlimits` was not given!";
		throw InputError( warning.str() );
	}

	// `Analysis` must be given if not given `--no-analysis`
	if ( !_given_xy && !_given_radial && !_no_analysis ){
		std::stringstream warning;
		warning << "In file `" << _rc_file << "`, the `Analysis` domain went ";
		warning << "unspecified, but the `--no-analysis` flag was not given!";
		throw InputError( warning.str() );
	}

	// check for conflicting assignments
	if ( _given_xy && _given_radial ){
		std::stringstream warning;
		warning << "In file `" << _rc_file << "`, there can be only one ";
		warning << "assignment for the `Analysis` domain!";
		throw InputError( warning.str() );
	}

	// set the `first_seed`
	convert.clear();
	convert.str( argument["--first-seed"] );
	if ( !given["--first-seed"] ) _first_seed = 19650218ULL;
	else if ( !(convert >> _first_seed) )
	throw InputError("--first-seed needs an interger value!");

	// check for `debug` mode
	_debug_mode = given["--debug"] ? true : false;
}

void Parser::SetLimits(const std::string &limits, const std::string &begin,
	const std::string &end){

	double numeric[2];
	std::stringstream convert(begin + " " + end);

	if ( !( convert >> numeric[0] ) || !( convert >> numeric[1] ) ){
		// these were not both numeric values!
		std::stringstream warning;
		warning << "In file `" << _rc_file << "` on line " << _line_number;
		warning << ", `" << limits << "` needs numeric values!";
		throw InputError( warning.str() );
	}

	std::vector<double> temp(numeric, numeric + 2);

	switch (limits[0]){

		case 'X': _x_limits = temp; _given_xlims = true; break;
		case 'Y': _y_limits = temp; _given_ylims = true; break;
		case 'Z': _z_limits = temp; _given_zlims = true; break;
		// uh oh, how did we get here?
		default: throw InputError("Parser::SetLimits() got non `XYZ` limit!");
	}
}

void Parser::SetXY(const std::string &xres, const std::string &yres){

	_analysis_domain = 3;
	double numeric[2];
	std::stringstream convert(xres + " " + yres);

	if ( !( convert >> numeric[0] ) || !( convert >> numeric[1] ) ){
		// these were not both numeric values!
		std::stringstream warning;
		warning << "In file `" << _rc_file << "` on line " << _line_number;
		warning << ", `XY` needs numeric values!";
		throw InputError( warning.str() );
	}

	std::vector<int> temp(numeric, numeric + 2);
	_xy_resolution = temp;
	_given_xy = true;
}

void Parser::SetRadial(const std::string &coord, const std::string &res){

	if (coord == "R")
		_analysis_domain = 1;
	else if (coord == "Rho")
		_analysis_domain = 2;
	else
		// how did we get here?
		throw InputError("Parser::SetRadial() got non R/Rho setting!");

	double numeric;
	std::stringstream convert(res);

	if ( !( convert >> numeric ) ){
		// this was not a numeric value!
		std::stringstream warning;
		warning << "In file `" << _rc_file << "` on line " << _line_number;
		warning << ", `" << coord << "` needs a numeric value!";
		throw InputError( warning.str() );
	}

	_radial_resolution = numeric;
	_given_radial = true;
}


void Parser::Set(const std::vector<std::string> &line){
	//
	// take a vector of cammands from the `_rc_file` and `Set` that
	// parameter.
	//

	// ensure that we have at least two more `words`
	if ( line.size() < 3 ){
		std::stringstream warning;
		warning << "In file `" << _rc_file << "` on line " << _line_number;
		warning << ", there are insufficient arguments!";
		throw InputError( warning.str() );
	}

	if ( argument.find(line[1]) != argument.end() ) {
		// user is assigning a runtime parameter
		argument[ line[1] ] = line[2];
		given[ line[1] ]    = true;
	}

	else if ( line[1] == "Xlimits" || line[1] == "Ylimits" ||
			line[1] == "Zlimits" ){
		// user is assigning limits for `the box`
		if ( line.size() < 4 ){
			// insufficient arguments
			std::stringstream warning;
			warning << "In file `" << _rc_file << "` on line " << _line_number;
			warning << ", `" << line[1] << "` requires two values!";
			throw InputError( warning.str() );
		}
		// assign the parameter
		SetLimits(line[1], line[2], line[3]);
	}

	else if ( line[1] == "Analysis" ){

		// check for the three allowed categories
		if ( line[2] != "R" && line[2] != "Rho" && line[2] != "XY"){
			std::stringstream warning;
			warning << "In file `" << _rc_file << "` on line " << _line_number;
			warning << ", `" << line[2] << "` is not a recognized ";
			warning << "option for `Analysis`!";
			throw InputError( warning.str() );
		}

		// we need one more parameter if it's XY
		if ( line[2] == "XY" ){

			// check for arguments
			if ( line.size() < 5 ){
				std::stringstream warning;
				warning << "In file `" << _rc_file << "` on line " << _line_number;
				warning << ", `XY` needs two values!";
				throw InputError( warning.str() );
			}

			SetXY(line[3], line[4]);
		}

		else {

			// check for arguments
			if ( line.size() < 4 ){
				std::stringstream warning;
				warning << "In file `" << _rc_file << "` on line " << _line_number;
				warning << ", `" << line[2] << "` needs a values!";
				throw InputError( warning.str() );
			}

			SetRadial(line[2], line[3]);
		}

	} else {
		// unrecognized command
		std::stringstream warning;
		warning << "In file `" << _rc_file << "` on line " << _line_number;
		warning << ", `" << line[1] << "` is not a recognized parameter!";
		throw InputError( warning.str() );
	}
}


void Parser::Include(const std::vector<std::string> &line){

	//
	// Parse a line of text from the RC-file for `include`ing a
	// profile
	//

	if ( line.size() < 2 ){

		std::stringstream warning;
		warning << "In file `" << _rc_file << "` on line " << _line_number;
		warning << ", you didn't specify a name of a profile!";
		throw InputError( warning.str() );
	}

	if ( line.size() > 3 ){

		// I don't know what your other argument was for!
		std::stringstream warning;
		warning << "In file `" << _rc_file << "` on line " << _line_number;
		warning << ", there are too many arguments! If you have spaces in ";
		warning << "your file path be sure to put it in quotes.";
		throw InputError( warning.str() );
	}

	if ( line.size() < 3){

		// we don't have the name of a file so we will just add the profile
		// without one.
		UsedPDFs[ line[1] ] = "";

	} else {

		UsedPDFs[ line[1] ] = line[2];
	}
}

// remove all characters after `delim`
void Parser::Clip(std::string &input_string, const std::string &delim){
    
    std::size_t pos = input_string.find(delim);
    
    if ( pos != std::string::npos )
        input_string.replace(pos, input_string.length(), "");
}

// split a string into `words`, grouping quoted words
std::vector<std::string> Parser::Split(const std::string &input){
    
    std::vector<std::string> sentence;
    std::string word;
    bool quoted = false;
    
    for (int i = 0; i < input.length(); i++){
        
        if (input[i] == ' '){
            
            if (quoted) word += " ";
            else {
                
                if ( !word.empty() )
                    sentence.push_back(word);
                
                word = "";
            }
            
        }
        
        else if ( input[i] == '\"' ) quoted = !quoted;
        else word += input[i];
    }
    
    // get last word
    if ( !word.empty() )
        sentence.push_back(word);
    
    return sentence;
}
    
// replace all instances of `search_str` in `input_str` with `replace_str`
void Parser::ReplaceAll(const std::string &search_str,
                const std::string &replace_str, std::string& input_str ){
    
    std::size_t pos = 0;
    
    while ( ( pos = input_str.find(search_str, pos) ) != std::string::npos ){
        
        input_str.replace( pos, search_str.length( ), replace_str );
        pos += replace_str.length( );
    }
}

// retrieval functions, `getters`
unsigned long long Parser::GetNumParticles() const {
	return _num_particles;
}

int Parser::GetNumTrials() const {
	return _num_trials;
}

int Parser::GetNumThreads() const {
	return _num_threads;
}

int Parser::GetVerbosity() const {
	return _verbose;
}

bool Parser::GetKeepPosFlag() const {
	return _keep_pos;
}

bool Parser::GetKeepRawFlag() const {
	return _keep_raw;
}

bool Parser::GetAnalysisFlag() const {
	return _no_analysis;
}

bool Parser::GetDebuggerFlag() const {
	return _debug_mode;
}

std::vector<double> Parser::GetXlimits() const {
	return _x_limits;
}

std::vector<double> Parser::GetYlimits() const {
	return _y_limits;
}

std::vector<double> Parser::GetZlimits() const {
	return _z_limits;
}

std::vector<int> Parser::GetXYResolution() const {
	return _xy_resolution;
}

int Parser::GetRadialResolution() const {
	return _radial_resolution;
}

int Parser::GetAnalysisDomain() const {
	return _analysis_domain;
}

std::string Parser::GetOutPath() const {
	return _out_path;
}

std::string Parser::GetRawPath() const {
	return _raw_path;
}

std::string Parser::GetTmpPath() const {
	return _tmp_path;
}

std::string Parser::GetPosPath() const {
	return _pos_path;
}

std::string Parser::GetRCFile() const {
	return _rc_file;
}

unsigned long long Parser::GetFirstSeed() const {
	return _first_seed;
}

std::map<std::string, std::string> Parser::GetUsedPDFs() const {
	return UsedPDFs;
}

} // namespace Gaia
