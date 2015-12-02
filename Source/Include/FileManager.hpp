// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// GNU General Public License v3.0
// Include/FileManager.hpp
//
// TODO: header

#ifndef _FILEMANAGER_HH_
#define _FILEMANAGER_HH_

#include <string>
#include <vector>
#include <map>

#include <Vector.hpp>

namespace Gaia {

class FileManager {

public:

    static FileManager* GetInstance();
    static void Release();
    ~FileManager(){}

    void Initialize();

    void SavePositions(const std::vector<Vector>&, const std::size_t );

    void SaveRaw( const std::vector<double>&, const std::size_t );

    void SaveMap( const std::map<std::string, std::vector<double>> Axis );

    void SaveOutput( const std::vector<double>&, const std::vector<double>&,
        const std::size_t );

    void SaveOutput( const std::vector< std::vector<double> >&,
        const std::vector< std::vector<double> >&, const std::size_t );

private:

    static FileManager* instance;
    FileManager(){}

    int verbose;
    std::string pos_path, raw_path, out_path, map_path;

};

} // namespace Gaia

#endif
