/*
* Filename : RemoveAmpliconDuplicates.h
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Status: Release
*/

#include <vector>
#include <unordered_map>
#include <string>

using namespace std;

//shared variable types
	typedef struct {
		string HeaderR1;
		string HeaderR2;
		string SeqR1;
		string SeqR2;
		string QualR1;
		string QualR2;
		double ReadErrors;
		double RTIErrors;
		unsigned long Frequency;
		bool PrintRead;
	} molecule;

	typedef struct {
		string HeaderR1;
		string HeaderR2;
		string SeqR1;
		string SeqR2;
		string QualR1;
		string QualR2;
	} unfilteredread;

	typedef struct {
		string AmpliconID;
		string FPrimer;
		string RPrimer;
		unsigned short FPrimerLen;
		unsigned short RPrimerLen;
	} amplicon;

	typedef struct {
		string RTI;
		bool Print;
		unsigned long Frequency;
	} usableRTI;

	//shared funtions
	double CalcReadErrorRate(const string& Qual, const unsigned QScorePhredOffset);
	bool getAmplicons(ifstream& AmpliconsIn, vector<amplicon>& Amplicons, unordered_map<string, bool>& AmpliconStrand);
	unsigned getHammingDistance(const string& str1, const string& str2);
	bool MatchPrimer(const string& Seq, const string& Primer);
	string ReverseComplement(const string& DNA);
	void RightPrimerClipper(string& Seq, string& Qual, const string& Primer);
	void FilterRTIsbyEditDistance(unordered_map<string, unordered_map<string, molecule>> & Reads, const unsigned MinRTIEditDistance);
	bool RTIQfilter(const string& Qual, const unsigned RTILen, const unsigned QScorePhredOffset, const unsigned MinRTIBaseQScore);
	string getSampleID(const string& FASTQFilename);
	void RTIDepthErrorRateFilter(unordered_map<string, unordered_map<string, molecule>> & Reads, const unsigned MinRTIDepthErrorRate);
	double getHighestErrorRate(const string& Qual, const unsigned QScorePhredOffset);

	void PrintParameters(int argc, char* argv[], const float ProgramVersion, const unsigned RTILen,
		const unsigned AntiComplementaryRegionLen, const unsigned MinRTIBaseQScore, const unsigned MinRTIEditDistance,
		const unsigned QScorePhredOffset, const unsigned MaxQScore, const unsigned MinInsertSize, const unsigned MinRTIDepthErrorRate);
	
	bool ReadMerger(const string& SeqR1, const string& QualR1, string SeqR2, string QualR2,
		const unsigned MaxQScore, const unsigned QScorePhredOffset, pair<string, string>& MergedRead);

	molecule MakeTempRead(const string& HeaderR1, const string& HeaderR2, const string& SeqR1, const string& SeqR2,
		const string& QualR1, const string& QualR2, const double& ReadErrors, const double& RTIErrors, const unsigned long& Frequency);
