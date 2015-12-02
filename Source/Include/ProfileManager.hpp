// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// GNU General Public License v3.0
// Include/ProfileManager.hpp
//

#ifndef _PROFILEMANAGER_HH_
#define _PROFILEMANAGER_HH_

#include <vector>

#include <ProfileBase.hpp>

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
