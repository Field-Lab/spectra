#ifndef __RAWDATAHEADER512__
#define __RAWDATAHEADER512__

#include <stdio>
#include <stdlib>
#include <string>
#include "RawDataFile.h"

class RawDataHeader512 {
public:
	static const bool DEBUG = false;

	static const int DEFAULT_TIME_BASE = 1904;
	
	static const int HEADER_LENGTH_TAG = 0;
	static const int TIME_TAG = 1;
	static const int COMMENT_TAG = 2;
	static const int FORMAT_TAG = 3;
	static const int ARRAY_ID_TAG = 4;
	static const int FREQUENCY_TAG = 5;
	static const int TRIGGER_TAG = 6;
	static const int DATASET_IDENTIFIER_TAG = 7;
	static const int DATA_TAG = 499;
	static const int FILE_TYPE = 0x512;
	
	RawDataHeader512 (long secondsTime, int nElectrodes, int frequency,
		int nSamples, int arrayID, int format, string datasetIdentifier,
		string comment);

	RawDataHeader512 ( int timeBase, long secondsTime, int nElectrodes,
		int frequency, int nSamples, int arrayID, int format,
		string datasetIdentifier, string comment );


	~RawDataHeader512 ();
};

#endif // __RAWDATAHEADER512__