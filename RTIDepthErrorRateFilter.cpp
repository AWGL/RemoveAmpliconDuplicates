/*
* Filename : RTIDepthErrorRateFilter.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Excludes low depth high error rate RTIs
* Status: Release
*/

#include <unordered_map>
#include <RemoveAmpliconDuplicates.h>

using namespace std;

void RTIDepthErrorRateFilter(unordered_map<string, unordered_map<string, molecule>> & Reads, const unsigned MinRTIDepthErrorRate){ //amplicon, RTI, molecule

	double AvgRTIErrorRate;

	//iterate over amplicons with usable reads
	for (auto & Amplicon : Reads){

		//iterate over RTIs associated with this amplicon
		for (auto & RTI : Amplicon.second){

			//skip over RTIs that will not be printed
			if (RTI.second.PrintRead == false ){
				continue;
			}

			AvgRTIErrorRate = (double) RTI.second.RTIErrors / RTI.second.Frequency;

			if (RTI.second.Frequency / AvgRTIErrorRate < MinRTIDepthErrorRate){
				RTI.second.PrintRead = false;
			}
		}
	}

	return;
}