// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// Include/ProfileManager.hh
//

#ifndef _PROFILEMANAGER_HH_
#define _PROFILEMANAGER_HH_

namespace GAIA {

class ProfileManager {

public:
	
	ProfileManager();
	~ProfileManager();
	
	Initialize();
	
	// maintain list of `known` and `used` PDFs
	ProfileBase *knownPDFs, *usedPDFs;
	
private:

};

} // namespace GAIA

#endif