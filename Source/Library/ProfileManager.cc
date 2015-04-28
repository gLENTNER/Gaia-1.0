// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Library/ProfileManager.cc
//
// TODO: source

#include <sstream>
#include <vector>
#include <string>
#include <map>

#include <ProfileManager.hh>
#include <Profiles.hh>

namespace Gaia {

ProfileManager::ProfileManager(){}

ProfileManager::~ProfileManager(){

    for ( auto& pdf : UsedPDFs ){

        ProfileBase *this_pdf = pdf;

        if( this_pdf ){
            delete this_pdf;
            this_pdf = nullptr;
        }
    }

    UsedPDFs.clear();
}

void ProfileManager::Initialize(){
    //
    // establish known profiles and declare what profiles will be
    // in use
    //

    // hard code known profiles
    std::vector<ProfileBase*> KnownPDFs;
    KnownPDFs.push_back( new MassDensity() );
    KnownPDFs.push_back( new Spiral() );
    KnownPDFs.push_back( new Metallicity() );
    KnownPDFs.push_back( new Habitability() );
    KnownPDFs.push_back( new Surface() );

    // map profiles to their names
    std::map<std::string, ProfileBase*> available;
    for ( auto& pdf : KnownPDFs )
        available[ pdf -> Name() ] = pdf;

    // grab the parser
    Parser *parser = Parser::GetInstance();

    // retrieve the map of used profiles from the parser
    std::map<std::string, std::string> given = parser -> GetUsedPDFs();

    for ( auto& profile : given ){

        if ( available.find( profile.first ) == available.end() ){

            std::stringstream warning;
            warning << "From file `" << parser -> GetRCFile() << "`, ";
            warning << "the given profile name `" << profile.first;
            warning << "` does not match any available. Be sure that the ";
            warning << "profiles defined in `Profiles.hh` match those ";
            warning << "listed in `ProfileManager::Initialize()` from ";
            warning << "`ProfileManager.cc`. \n";
            throw ProfileError( warning.str() );
        }

        if ( profile.second.empty() ){

            // this is an analytical profile (no need for initialization)
            UsedPDFs.push_back( available[profile.first] );

        } else {

            // we must initialize it first
            ProfileBase *this_pdf = available[profile.first];
            this_pdf -> Initialize( profile.second );
            UsedPDFs.push_back( this_pdf );
        }

    }
}

} // namespace Gaia
