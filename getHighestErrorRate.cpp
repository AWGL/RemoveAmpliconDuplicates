/*
* Filename : getHighestErrorRate.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Returns highest probability of base error across entire sequence
* Status: Release
*/

#include <string>
#include <math.h>
#include <RemoveAmpliconDuplicates.h>

using namespace std;

double getHighestErrorRate(const string& Qual, const unsigned QScorePhredOffset){

	double HighestBaseError = 0;
	unsigned phredScore;

	for (unsigned n = 0; n < Qual.length(); ++n){

		//get PhredScore 1-40
		phredScore = Qual[n] - QScorePhredOffset;

		if ((double)pow(10.00, (double)phredScore / -10.00) > HighestBaseError){
			HighestBaseError = (double)pow(10.00, (double)phredScore / -10.00);
		}

	}

	return HighestBaseError;

}