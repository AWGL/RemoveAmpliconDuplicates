/*
* Filename : CalcReadErrorRate.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Returns average probability of base error across entire read
* Status: Release
*/

#include <string>
#include <math.h>
#include <RemoveAmpliconDuplicates.h>

using namespace std;

double CalcReadErrorRate(const string& Qual, const unsigned QScorePhredOffset){ //implementation of: http://www.drive5.com/usearch/manual/avgq.html

	double errorRate = 0;
	unsigned phredScore;

	for (unsigned n = 0; n < Qual.length(); ++n){

		//get PhredScore 1-40
		phredScore = Qual[n] - QScorePhredOffset;

		errorRate += (double)pow(10.00, (double)phredScore / -10.00);

	}

	return errorRate / Qual.length();

}