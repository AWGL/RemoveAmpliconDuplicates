/*
* Filename : MatchPrimer.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Uses Smith-Waterman (SeqAn) local alignement to identify supplied primer sequences within the read
* Status: Release
*/

#include <string>
#include <seqan/align.h>
#include <RemoveAmpliconDuplicates.h>

using namespace std;

bool MatchPrimer(const string& Seq, const string& Primer) //iterate over bases of primer and match to seq
{
	seqan::Align< seqan::String<char> > Alignment;
	seqan::resize(rows(Alignment), 2); //pairwise
	seqan::assignSource(row(Alignment, 0), Seq);
	seqan::assignSource(row(Alignment, 1), Primer);

	//match mismatch extend open
	if (seqan::localAlignment(Alignment, seqan::Score<int>(1, -2, -4)) >= 10 &&
		seqan::clippedBeginPosition(row(Alignment, 0)) == 0 &&
		seqan::clippedBeginPosition(row(Alignment, 1)) == 0){ //mapQ & left most position

		return 1;

	} else {
		return 0;
	}

}

/*bool MatchPrimer(string& seq, string& primer) //iterate over bases of primer and match to seq
{
	float BasesMatched = 0;
	unsigned MaxMismatchLen = 3, PrimerLen = primer.length(); //no mismatches in the last 3bp -- prevents indels through phase shift and reduced off-target reads

	for (unsigned base = 0; base < PrimerLen; ++base) {

		if (seq[base] == primer[base]) {
			BasesMatched++;
		} else if (base > (PrimerLen - MaxMismatchLen)) {
			return 0;
		}

	}

	if (BasesMatched / (PrimerLen - MaxMismatchLen) > 0.8) { //check if match is acceptable
		return 1;
	}

	return 0;
}*/