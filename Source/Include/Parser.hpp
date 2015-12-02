// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// GNU General Public License v3.0
// Include/Parser.hpp
//
// This header file contains the class template for the `Parser` object.
// The purpose of this singleton class is to interpret the arguments passed
// to `main` at runtime and to read the configuration file, retaining the
// parameters for retrieval by the other objects.

#ifndef _PARSER_HH_
#define _PARSER_HH_

#include <vector>
#include <string>
#include <map>

namespace Gaia {

class Parser {

public:

	static Parser* GetInstance();
	static void Release();
	~Parser(){}

	// set up the simulation environment
	void Setup(const int argc, const char *argv[]);

	// retrieval functions, "getters"
	std::size_t GetNumParticles() const;
	int GetNumTrials() const;
	int GetNumThreads() const;
	int GetVerbosity() const;
	bool GetKeepPosFlag() const;
	bool GetKeepRawFlag() const;
	bool GetAnalysisFlag() const;
	bool GetDebuggerFlag() const;
	std::vector<double> GetXlimits() const;
	std::vector<double> GetYlimits() const;
	std::vector<double> GetZlimits() const;
	std::vector<std::size_t> GetResolution() const;
	std::vector<std::string> GetAxes() const;
	std::string GetOutPath() const;
	std::string GetRawPath() const;
	std::string GetPosPath() const;
	std::string GetMapPath() const;
	std::string GetRCFile() const;
	unsigned long long GetFirstSeed() const;
	double GetSampleRate() const;
	double GetMeanBandwidth() const;
	double GetStdevBandwidth() const;
	std::map<std::string, std::string> GetUsedPDFs() const;

private:

	static Parser* instance;
	Parser(){}

	// set default values for `argument` map, etc...
	void SetDefaults();

	// parse the commands in the configuration (rc) file
	void ReadRC();

	// parse the input arguments from `main`
	void Interpret();

	// check/assign arguments
	void Rectify();

	// helper functions
	void SetLimits(const std::string&, const std::string&, const std::string&);
	void SetAnalysis(const std::vector<std::string>& );

	void Clip(std::string &input_string, const std::string &delim);
	std::vector<std::string> Split(const std::string &input);
	void ReplaceAll(const std::string &search_str,
	        const std::string &replace_str, std::string& input_str);

	// helper functions for `Command` map
	void Set(const std::vector<std::string>&);
	void Include(const std::vector<std::string>&);

	// maps of parameters
	std::map<std::string, std::string> argument, implicit;
	std::map<std::string, bool> given;

	// simulation parameters, see SetDefaults() for defaults
	int _verbose, _num_threads, _num_trials, _line_number;
	bool _keep_raw, _keep_pos, _no_analysis, _debug_mode;
	std::size_t _num_particles;
	std::string _out_path, _raw_path, _pos_path, _map_path, _rc_file;
	unsigned long long _first_seed;
	double _sample_rate, _mean_bandwidth, _stdev_bandwidth;

	// vector for `argv`
	std::vector<std::string> _cmd_args;

	// items from rc file
	std::vector<double> _x_limits, _y_limits, _z_limits;
	std::vector<std::size_t> _resolution;
	std::vector<std::string> _axes;

	bool _given_xlims, _given_ylims, _given_zlims, _given_analysis;

	// map of profile names from RC-file with file paths
	std::map<std::string, std::string> UsedPDFs;
};

} // namespace Gaia

#endif
