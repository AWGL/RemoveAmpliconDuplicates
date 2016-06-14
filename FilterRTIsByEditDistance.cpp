/*
* Filename : FilterRTIsbyEditDistance.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Sorts supplied vector high to low and marks random template identifiers for discard if they have an edit distance less than specified
* Status: Release
*/

#include <unordered_map>
#include <RemoveAmpliconDuplicates.h>

using namespace std;

void FilterRTIsbyEditDistance(unordered_map<string, unordered_map<string, molecule>> & Reads, const unsigned MinRTIEditDistance){ //amplicon, RTI, molecule

	unsigned HammingDistance;

	//iterate over amplicons with usable reads
	for (auto & Amplicon : Reads){

		//iterate over RTIs associated with this amplicon
		for (auto & OuterRead : Amplicon.second){

			//skip over RTIs that will not be printed
			if (OuterRead.second.PrintRead == false){
				continue;
			}

			//iterate back over over RTIs associated with this amplicon
			for (auto & InnerRead : Amplicon.second){

				//skip over RTIs that will not be printed
				if (InnerRead.second.PrintRead == false){
					continue;
				}

				//calculate edit distance
				HammingDistance = getHammingDistance(OuterRead.first, InnerRead.first);

				if (HammingDistance != 0 && HammingDistance < MinRTIEditDistance){ //too similar discard RTI

					//retain highest frequency RTI
					if (OuterRead.second.Frequency > InnerRead.second.Frequency){
						InnerRead.second.PrintRead = false; // this record will not be printed
					} else {
						OuterRead.second.PrintRead = false; // this record will not be printed
					}

				}

			}

		}
	}

	return;
}