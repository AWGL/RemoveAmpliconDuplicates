/*
* Filename : MakeTempRead.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Returns TempRead structure
* Status: Release
*/

#include <string>
#include <RemoveAmpliconDuplicates.h>

using namespace std;

molecule MakeTempRead(const string& HeaderR1, const string& HeaderR2, const string& SeqR1, const string& SeqR2,
	const string& QualR1, const string& QualR2, const double& ReadErrors, const double& RTIErrors, const unsigned long& Frequency) {

	molecule TempRead;

	TempRead.ReadErrors = ReadErrors;
	TempRead.RTIErrors = RTIErrors;
	TempRead.Frequency = Frequency;
	TempRead.HeaderR1 = HeaderR1;
	TempRead.HeaderR2 = HeaderR2;
	TempRead.QualR1 = QualR1;
	TempRead.QualR2 = QualR2;
	TempRead.SeqR1 = SeqR1;
	TempRead.SeqR2 = SeqR2;
	TempRead.PrintRead = true;

	return TempRead;
}