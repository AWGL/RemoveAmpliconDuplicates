/*
* Filename : PrintParameters.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Prints program parameters to user for logging
* Status: Release
*/

#include <iostream>
#include <RemoveAmpliconDuplicates.h>

using namespace std;

void PrintParameters(int argc, char* argv[], const float ProgramVersion, const unsigned RTILen, 
	const unsigned AntiComplementaryRegionLen, const unsigned MinRTIBaseQScore, const unsigned MinRTIEditDistance,
	const unsigned QScorePhredOffset, const unsigned MaxQScore, const unsigned MinInsertSize, const unsigned MinRTIDepthErrorRate){

	//print parameters
	cout << "\nRemoveAmpliconDuplicates v" << ProgramVersion << endl;

	//print commandline arguments
	cout << "CL:";
	for (unsigned n = 0; n < argc; ++n){
		cout << ' ' << argv[n];
	}

	cout << "\n" << endl;

	//print varaibles
	cout << "RTILength: " << RTILen << endl;
	cout << "AntiComplementaryRegionLength: " << AntiComplementaryRegionLen << endl;
	cout << "MinimumRTIBaseQScore: " << MinRTIBaseQScore << endl;
	cout << "MinimumRTIEditDistance: " << MinRTIEditDistance << endl;
	cout << "MinRTIDepthErrorRate: " << MinRTIDepthErrorRate << endl;
	cout << "QScorePhredOffset: " << QScorePhredOffset << endl;
	cout << "MaxQScore: " << MaxQScore << endl;
	cout << "MinInsertSize: " << MinInsertSize << endl;

	return;
}