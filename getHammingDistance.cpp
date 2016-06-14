/*
* Filename : getHammingDistance.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Returns hamming distance of two strings of the same length.
* Status: Release
*/

#include <string>
#include <RemoveAmpliconDuplicates.h>

using namespace std;

unsigned getHammingDistance(const string& str1, const string& str2){

	unsigned HammingDistance = 0;

	for (unsigned n = 0; n < str1.length(); ++n){ //must be the same length
		if (str1[n] != str2[n]){
			HammingDistance++;
		}
	}

	return HammingDistance;
}
