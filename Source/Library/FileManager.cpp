// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Library/FileManager.cc
//
// #TODO:0 source

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <FileManager.hpp>
#include <Parser.hpp>
#include <Vector.hpp>
#include <Exception.hpp>

namespace Gaia {

// static pointer for singleton class
FileManager* FileManager::instance = nullptr;

// retrieval of pointer
FileManager* FileManager::GetInstance() {

    if ( !instance )
        instance = new FileManager();

    return instance;
}

// destroy the singleton
void FileManager::Release( ) {

    if ( instance ) {
        delete instance;
        instance = nullptr;
    }
}

void FileManager::Initialize(){

    // get info from parser
    Parser *parser = Parser::GetInstance();
    verbose  = parser -> GetVerbosity();
    pos_path = parser -> GetPosPath();
    raw_path = parser -> GetRawPath();
    out_path = parser -> GetOutPath();
    map_path = parser -> GetMapPath();
}

void FileManager::SavePositions(const std::vector<Vector> &positions,
    const std::size_t trial ){

    // build file name
    std::stringstream buffer;
    buffer << pos_path << trial << ".dat";
    std::string filename = buffer.str();

    // open file and write positions
    std::ofstream output( filename.c_str() );

    if (output) {

        output.precision(16);

        if (verbose) std::cout
            << "\n\n Saving position vectors to `"
            << filename << "` ... ";
            std::cout.flush();

        for ( const auto& vec : positions )
            output << vec << std::endl;

        if (verbose)
            std::cout << "done\n";
            std::cout.flush();

    } else throw IOError("From FileManager::SavePositions(), I "
    "couldn't open the file `" + filename + "`!");
}

void FileManager::SaveRaw(const std::vector<double> &seperations,
    const std::size_t trial){

    // build file name
    std::stringstream buffer;
    buffer << raw_path << trial << ".dat";
    std::string filename = buffer.str();

    // open file and write positions
    std::ofstream output( filename.c_str() );

    if (output) {

        output.precision(16);

        if (verbose) std::cout
            << "\n\n Saving raw nearest neighbor distances to `"
            << filename << "` ... ";
            std::cout.flush();

        for ( const auto& sep : seperations )
            output << sep << std::endl;

        if (verbose)
            std::cout << "done\n";
            std::cout.flush();

    } else throw IOError("From FileManager::SaveRaw(), I "
    "couldn't open the file `" + filename + "`!");
}

void FileManager::SaveMap( const std::map<std::string, std::vector<double>> Axis ){

    //
    // Save the axis information the results were mapped to.
    //

    for ( const auto& ax : Axis ){

        // build file name
        std::stringstream buffer;
        buffer << map_path << ax.first << ".dat";
        std::string filename = buffer.str();

        if (verbose) std::cout
            << "\n Saving analysis coordinate map data to `"
            << filename << "` ... ";
            std::cout.flush();

        // open file and write out `map` coordinates
        std::ofstream output( filename.c_str() );

        if (output) {

            output.precision(16);

            for ( const auto& x : ax.second )
                output << x << std::endl;

        } else throw IOError("From FileManager::SaveMap(), I couldn't open "
        "the file, `" + filename + "`!");

        if (verbose)
            std::cout<< "done\n";
            std::cout.flush();
    }
}

void FileManager::SaveOutput(const std::vector<double> &mean,
    const std::vector<double> &variance, const std::size_t trial){

    //
    // Save mean and standard deviation to file
    //

    // build file name
    std::stringstream buffer;
    buffer << out_path << trial << ".dat";
    std::string filename = buffer.str();

    if (verbose) std::cout
        << "\n\n Saving Profile data to `"
        << filename << "` ... ";
        std::cout.flush();

    // take the sqrt of variance for the standard deviation
    std::vector<double> stdev(variance.size(), 0.0);
    for (std::size_t i = 0; i < stdev.size(); i++)
        stdev[i] = std::sqrt(variance[i]);

    // open file and write positions
    std::ofstream output( filename.c_str() );

    if (output) {

        output.precision(16);

        for (std::size_t i = 0; i < mean.size(); i++)
            output << mean[i] << " " << stdev[i] << std::endl;

    } else throw IOError("From FileManager::SaveOutput(), I "
        "couldn't open the file, `" + filename + "`!");

    if (verbose)
        std::cout << "done\n";
        std::cout.flush();
}

void FileManager::SaveOutput( const std::vector< std::vector<double> > &mean,
    const std::vector< std::vector<double> > &variance, const std::size_t trial){

        //
        // Save mean and standard deviation to file, 2D
        //

        // take the sqrt of variance for the standard deviation
        std::vector< std::vector<double> > stdev = variance;
        for ( auto& row : stdev )
        for ( auto& element : row )
                element = std::sqrt(element);

        // create map to shorten notation ...
        std::map< std::string, std::vector< std::vector<double> > > results;
        results["mean"]  = mean;
        results["stdev"] = stdev;

        // counting variable updates filename
        for( const auto& matrix : results ){

                // build file name
                std::stringstream buffer;
                buffer << out_path << trial << "-" << matrix.first << ".dat";
                std::string filename = buffer.str();

                if (verbose) std::cout
                        << "\n\n Saving `" << matrix.first << "` data to file, `"
                        << filename << "` ... ";
                        std::cout.flush();

                // open file and write positions
                std::ofstream output( filename.c_str() );

                if (output) {

                        output.precision(16);

                        for ( const auto& row : matrix.second ){

                                for ( const auto& element : row )
                                        output << element << " ";

                                output << std::endl;
                        }

                } else throw IOError("From FileManager::SaveOutput(), I "
                        "couldn't open the file, `" + filename + "`!");

                if (verbose)
                        std::cout << "done";
                        std::cout.flush();
        }
}

} // namespace Gaia
