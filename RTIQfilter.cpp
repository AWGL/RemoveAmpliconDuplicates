/*
* Filename : RTIQfilter.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Returns true if all bases of the random template identifier are above the specified quality
* Status: Release
*/

#include <string>
#include <RemoveAmpliconDuplicates.h>

using namespace std;

bool RTIQfilter(const string& Qual, const unsigned RTILen, const unsigned QScorePhredOffset, const unsigned MinRTIBaseQScore){ //Check all bases of RTI are above minQx

	for (unsigned n = 0; n < RTILen; ++n){

		if (Qual[n] - QScorePhredOffset < MinRTIBaseQScore){
			return false;
		}

	}

	return true;
}