/*
* Filename : getAmplicons.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Extracts information from the supplied amplicon input file
* Status: Release
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <unordered_map>
#include <RemoveAmpliconDuplicates.h>

using namespace std;

bool getAmplicons(ifstream& AmpliconsIn, vector<amplicon>& Amplicons, unordered_map<string, bool>& AmpliconStrand){ //return success or failure

	unsigned short n, j = 0;
	string AmpliconLine;
	amplicon Temp;
	vector<string> AmpliconFields;

	if (AmpliconsIn.is_open()) {
		while (AmpliconsIn.good()) {
			getline(AmpliconsIn, AmpliconLine);

			AmpliconFields.clear();

			boost::trim(AmpliconLine); //remove whitespace at either end of line; MS excel likes putting this in.

			//skip empty lines and headers
			if (AmpliconLine == "" || AmpliconLine[0] == '#') {
				continue;
			} else {

				//tokenize string
				boost::split(AmpliconFields, AmpliconLine, boost::is_any_of("\t"), boost::token_compress_on); //substrings in elements

				if (AmpliconFields.size() != 4) {

					cerr << "ERROR: Amplicon list improperly formatted." << endl;
					cerr << "AmpliconList: AmpliconID ForwardPrimer ReversePrimer\n" << endl;
					return 1;

				} else {

					boost::to_upper(AmpliconFields[1]); //convert sequence to upper-case
					boost::to_upper(AmpliconFields[2]); //convert sequence to upper-case

					//check forward primer is standard DNA seq
					for (n = 0; n < AmpliconFields[1].length(); ++n){
						if (AmpliconFields[1][n] != 'A' &&
							AmpliconFields[1][n] != 'T' &&
							AmpliconFields[1][n] != 'G' &&
							AmpliconFields[1][n] != 'C'){

							cerr << "ERROR: Forward primer contains non-standard bases (only A,T,G or C allowed)" << endl;
							return 1;

						}
					}

					//check reverse primer is standard DNA seq
					for (n = 0; n < AmpliconFields[2].length(); ++n){
						if (AmpliconFields[2][n] != 'A' &&
							AmpliconFields[2][n] != 'T' &&
							AmpliconFields[2][n] != 'G' &&
							AmpliconFields[2][n] != 'C'){

							cerr << "ERROR: Reverse primer contains non-standard bases (only A,T,G or C allowed)" << endl;
							return 1;

						}
					}

					//bank primers
					Temp.AmpliconID = AmpliconFields[0];
					Temp.FPrimer = AmpliconFields[1];
					Temp.RPrimer = AmpliconFields[2];
					Temp.FPrimerLen = AmpliconFields[1].length();
					Temp.RPrimerLen = AmpliconFields[2].length();
					
					//get amplicon strand
					if (AmpliconFields[3] == "1"){
						AmpliconStrand[Temp.AmpliconID] = 1;
					} else if (AmpliconFields[3] == "0"){
						AmpliconStrand[Temp.AmpliconID] = 0;
					} else {
						cerr << "ERROR: Strand field must contain 0 or 1." << endl;
						return 1;
					}

					Amplicons.push_back(Temp);

				}

			}

		}

		AmpliconsIn.close();

	} else {
		cerr << "ERROR: Unable to open amplicon file" << endl;
		return 1;
	}

	return 0;

}