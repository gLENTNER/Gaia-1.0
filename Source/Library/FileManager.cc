// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Library/FileManager.cc
//
// TODO: source

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <FileManager.hh>
#include <Parser.hh>
#include <Vector.hh>
#include <Exception.hh>

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
    tmp_path = parser -> GetTmpPath();
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
        
        for ( const auto& vec : positions )
            output << vec << std::endl;
        
    } else throw IOError("From FileManager::SavePositions(), I "
    "couldn't open the file `" + filename + "`!");
}

} // namespace Gaia


































