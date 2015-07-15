#ifndef __RAWDATAFILE__
#define __RAWDATAFILE__

#include <iostream>
#include <fstream>
#include "RawDataHeader512.h"

class RawDataFile {

private:
	RawDataHeader512 header;

	RawDataFile ();
	
	~RawDataFile ();
};

#endif // __RAWDATAFILE__

