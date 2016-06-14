#include <string>
#include <RemoveAmpliconDuplicates.h>

using namespace std;

string getSampleID(const string& FASTQFilename){

	string SampleID;
	bool PassedFirstUnderScore = false;

	//iterate over FASTQFilename
	for (unsigned n = 0; n < FASTQFilename.size(); ++n){

		if (FASTQFilename[n] != '_'){
			SampleID += FASTQFilename[n];
		} else if (PassedFirstUnderScore == false){ //on first underscore
			SampleID += FASTQFilename[n];
			PassedFirstUnderScore = true;
		} else if (PassedFirstUnderScore == true){ //on second underscore
			break;
		}

	}

	return SampleID;

}