// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/FileManager.hh
//
// TODO: header

#ifndef _FILEMANAGER_HH_
#define _FILEMANAGER_HH_

#include <string>
#include <vector>

#include <Vector.hh>

namespace Gaia {
    
class FileManager {

public:
    
    static FileManager* GetInstance();
    static void Release();
    ~FileManager(){}
    
    void Initialize();
    void SavePositions(const std::vector<Vector>&, const std::size_t );
    void SaveRaw( const std::vector<double>&, const std::size_t );
    
private:
    
    static FileManager* instance;
    FileManager(){}
    
    int verbose;
    std::string pos_path, raw_path, out_path, tmp_path;
    
};
    
} // namespace Gaia

#endif