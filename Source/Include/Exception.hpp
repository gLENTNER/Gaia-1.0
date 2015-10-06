// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/Exception.hpp
//
// This header file contains a base and derived exception classes
// specific to the GAIA project.

#ifndef _EXCEPTION_HH_
#define _EXCEPTION_HH_

#include <exception>
#include <string>

namespace Gaia {

// base exception class for GAIA project
class Exception : public std::exception {

public:

	explicit Exception(const std::string& msg): _msg(msg){ }
	virtual ~Exception() throw() { }
	virtual const char* what() const throw(){ return _msg.c_str(); }

protected:

	std::string _msg;
};

// exception related to input/output file writing
class IOError : public Exception {
public:
	IOError( const std::string& msg ): Exception(
		"\n --> IOError: " + msg) { }
};

// exception from GAIA::Parser()
class InputError : public Exception {
public:
	InputError( const std::string& msg ): Exception(
		"\n --> InputError: " + msg) { }
};

// exception specific to the construction of `Profile`s
class ProfileError : public Exception {
public:
	ProfileError( const std::string& msg ): Exception(
		"\n --> ProfileError: " + msg){ }
};

// DivError is for Vector::Phi(), I can't assume what the polar
// angle is with no magnitude?
class DivError : public Exception {
public:
	DivError(const std::string& msg): Exception(
		"\n --> DivError: " + msg){ }
};

// index out of bounds error
class IndexError : public Exception {
public:
	IndexError(const std::string& msg): Exception(
		"\n --> IndexError: " + msg){ }
};

// not an error, but halt execution
class Usage : public Exception {
public:
	Usage(const std::string& msg): Exception("\n usage: " + msg){ }
};

} // namespace Gaia

#endif
