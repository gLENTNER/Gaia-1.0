// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/ProfileManager.hh
//

#ifndef _PROFILEMANAGER_HH_
#define _PROFILEMANAGER_HH_

#include <vector>

#include <ProfileBase.hh>

namespace Gaia {

class ProfileManager {

public:
	
    ProfileManager();
	~ProfileManager();
	
	void Initialize();
    
    // maintain list `used` PDFs
    std::vector<ProfileBase*> UsedPDFs;

};

} // namespace Gaia

#endif