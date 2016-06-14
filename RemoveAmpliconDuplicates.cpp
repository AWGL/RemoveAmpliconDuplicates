/*
* Filename : RemoveAmpliconDuplicates.cpp
* Author : Matthew Lyon, Wessex Regional Genetics Laboratory, Salisbury, UK & University of Southampton, UK
* Contact : mlyon@live.co.uk
* Description : Removes PCR duplicates from paired-end sequenced amplicon library preps using dual-tagged random template identifiers
* Status: Release
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <RemoveAmpliconDuplicates.h>

using namespace std;

int RandomGenerator(int i) {

	random_device rnd;
	return rnd() % i;

}

int main(int argc, char* argv[]) {

	const float ProgramVersion = 0.4;

	//check argument number is correct; print usage
	if (argc != 4) { //program ampliconlist r1 r2
		cerr << "\nProgram: RemoveAmpliconDuplicates v" << ProgramVersion << ' ' << __DATE__ << ' ' << __TIME__ << endl;
		cerr << "Contact: Matthew Lyon, WRGL/UoS (mlyon@live.co.uk)\n" << endl;
		cerr << "Usage: RemoveAmpliconDuplicates <AmpliconList> <R1.fastq> <R2.fastq>\n" << endl;
		cerr << "AmpliconList: AmpliconID ForwardPrimer ReversePrimer Strand\n" << endl;
		return -1;
	}

	//settings
	const unsigned RTILen = 5; //Random template identifier
	const unsigned AntiComplementaryRegionLen = 3;
	const unsigned MinRTIBaseQScore = 17; //Every base of the random template identifier (20% of counters contain one error)
	const unsigned MinRTIEditDistance = 2; //Minimum random template identifier edit distance
	const unsigned QScorePhredOffset = 33, MaxQScore = 40; //ILMN 1.8/1.9
	const unsigned MinInsertSize = 5;
	const unsigned MinRTIDepthErrorRate = 1000; //(RTI depth / highest base error rate of RTI) filter

	//stats
	unsigned long TotalPairedReads = 0, LenDiscardedReads = 0, RTIQualityDiscardedReads = 0,
		PrimerMatchedReads = 0, NMaskedReads = 0, TotalUsableMolecules = 0, TotalUsableReads = 0;
	unordered_map<string, unsigned long> AmpliconUsableReads; //total reads passing filter per amplicon
	unordered_map<string, unsigned long> AmpliconUniqueReads; //total reads after removing dups per amplicon

	//variables
	unsigned LineNo = 0, GetHeader = 0, HammingDistance, n;
	double ReadErrors, RTIErrors;
	molecule SavedBestRead, TempRead;
	unfilteredread TempUnfilteredRead;
	unordered_map<string, vector<unfilteredread>> UnfilteredReads;
	string FQLineR1, FQLineR2, HeaderR1, HeaderR2, SeqR1, SeqR2, QualR1, QualR2, RTI, RTIQualities;
	vector<amplicon> Amplicons;
	unordered_map<string, unordered_map<string, molecule>> Reads; //<ampliconID><RTI> = molecule
	unordered_map<string, bool> AmpliconStrand;

	//define input filenames & SampleID
	string AmpliconfN = argv[1], R1fN = argv[2], R2fN = argv[3];
	string SampleID = getSampleID(R1fN);

	//Open files for reading or writing
	ifstream AmpliconsIn(AmpliconfN.c_str());
	ifstream R1FQIn(R1fN.c_str());
	ifstream R2FQIn(R2fN.c_str());
	ofstream R1Dedupped0((R1fN + ".Dedupped_0.fastq").c_str(), ios::binary);
	ofstream R2Dedupped0((R2fN + ".Dedupped_0.fastq").c_str(), ios::binary);
	ofstream R1Trimmed0((R1fN + ".Trimmed_0.fastq").c_str(), ios::binary);
	ofstream R2Trimmed0((R2fN + ".Trimmed_0.fastq").c_str(), ios::binary);
	ofstream R1Dedupped1((R1fN + ".Dedupped_1.fastq").c_str(), ios::binary);
	ofstream R2Dedupped1((R2fN + ".Dedupped_1.fastq").c_str(), ios::binary);
	ofstream R1Trimmed1((R1fN + ".Trimmed_1.fastq").c_str(), ios::binary);
	ofstream R2Trimmed1((R2fN + ".Trimmed_1.fastq").c_str(), ios::binary);
	ofstream StatsOut((R1fN.substr(0, R1fN.find_first_of('_')) + "_RTIs.txt").c_str());
	ofstream RTIHeadersOut((R1fN.substr(0, R1fN.find_first_of('_')) + "_RTIHeaders.txt").c_str());

	//print stats headers
	StatsOut << "SampleID\tAmplicon\tStrand\tRTI\tFrequency (Reads)\tSequenceErrors\n";

	//print input pararmeters to user for logging
	PrintParameters(argc, argv, ProgramVersion, RTILen, AntiComplementaryRegionLen, MinRTIBaseQScore, 
		MinRTIEditDistance, QScorePhredOffset, MaxQScore, MinInsertSize, MinRTIDepthErrorRate);

	//store amplicon fields
	if (getAmplicons(AmpliconsIn, Amplicons, AmpliconStrand) == 1){
		return -1; //error with amplicon input
	}

	//parse FASTQs
	if (R1FQIn.is_open() && R2FQIn.is_open()) {
		while (R1FQIn.good() && R2FQIn.good()) {
		
			getline(R1FQIn, FQLineR1);
			getline(R2FQIn, FQLineR2);

			if (FQLineR1 == "" || FQLineR2 == "") { //Skip empty lines
				continue;
			}

			LineNo++;

			if (LineNo == 1) {

				HeaderR1 = FQLineR1;
				HeaderR2 = FQLineR2;

				if (GetHeader < 10){ //Check header hamming distance equals 1

					if (HeaderR1.length() != HeaderR2.length()){
						cerr << "ERROR: Read header hamming distance does not equal one. Check FASTQ input." << endl;
						return -1;
					} else if (getHammingDistance(HeaderR1, HeaderR2) != 1){
						cerr << "ERROR: Read header hamming distance does not equal one. Check FASTQ input." << endl;
						return -1;
					}

					GetHeader++;
				}

				TotalPairedReads++;

			} else if (LineNo == 2) {
				SeqR1 = FQLineR1;
				SeqR2 = FQLineR2;
			} else if (LineNo == 4) {
				QualR1 = FQLineR1;
				QualR2 = FQLineR2;
				
				LineNo = 0;

				if (SeqR1 == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" || SeqR2 == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"){
					NMaskedReads++;
					continue;
				}

				//Check if RTI consists of Qx bases
				if (RTIQfilter(QualR1, RTILen, QScorePhredOffset, MinRTIBaseQScore) == false || RTIQfilter(QualR2, RTILen, QScorePhredOffset, MinRTIBaseQScore) == false){
					RTIQualityDiscardedReads++;
					continue; //skip counters with any bases less than minQscore
				}

				//Define RTI
				RTI = SeqR1.substr(0, RTILen) + ReverseComplement(SeqR2.substr(0, RTILen));
				RTIQualities = QualR1.substr(0, RTILen) + QualR2.substr(0, RTILen);
					
				//Trim RTI
				SeqR1 = SeqR1.substr(RTILen + AntiComplementaryRegionLen, string::npos);
				SeqR2 = SeqR2.substr(RTILen + AntiComplementaryRegionLen, string::npos);
				QualR1 = QualR1.substr(RTILen + AntiComplementaryRegionLen, string::npos);
				QualR2 = QualR2.substr(RTILen + AntiComplementaryRegionLen, string::npos);

				//Iterate over amplicons
				for (n = 0; n < Amplicons.size(); ++n) {

					if (MatchPrimer(SeqR1, Amplicons[n].FPrimer) == 1) { //local alignment
						if (MatchPrimer(SeqR2, Amplicons[n].RPrimer) == 1) { //read matches to this amplicon

							PrimerMatchedReads++;

							//Trim adapter and right RTI from read
							RightPrimerClipper(SeqR1, QualR1, ReverseComplement(Amplicons[n].RPrimer));
							RightPrimerClipper(SeqR2, QualR2, ReverseComplement(Amplicons[n].FPrimer));

							//Reduce primer dimer; insert size less than MinInsertLength ignored
							if (SeqR1.length() > Amplicons[n].FPrimerLen + Amplicons[n].RPrimerLen + MinInsertSize &&
								SeqR2.length() > Amplicons[n].FPrimerLen + Amplicons[n].RPrimerLen + MinInsertSize) {

								AmpliconUsableReads[Amplicons[n].AmpliconID]++;
								TotalUsableReads++;

								//Add read pair to vector for downsampling
								TempUnfilteredRead.HeaderR1 = HeaderR1;
								TempUnfilteredRead.HeaderR2 = HeaderR2;
								TempUnfilteredRead.QualR1 = QualR1;
								TempUnfilteredRead.QualR2 = QualR2;
								TempUnfilteredRead.SeqR1 = SeqR1;
								TempUnfilteredRead.SeqR2 = SeqR2;

								UnfilteredReads[Amplicons[n].AmpliconID].push_back(TempUnfilteredRead);

								//print read headers associated with each RTI
								RTIHeadersOut << HeaderR1 << "\t" << RTI << "\n";

								//Calculate number of readErrors across both reads & RTI
								ReadErrors = CalcReadErrorRate(QualR1, QScorePhredOffset) + CalcReadErrorRate(QualR2, QScorePhredOffset);
								RTIErrors = getHighestErrorRate(RTIQualities, QScorePhredOffset);

								//check if this RTI has been seen before
								if (Reads[Amplicons[n].AmpliconID].count(RTI) == 1){ //amplicon:RTI = molecule

									SavedBestRead = Reads[Amplicons[n].AmpliconID][RTI]; //get current best read for this RTI

									if (SavedBestRead.ReadErrors > ReadErrors){ //overwrite old read with new read containing less readErrors

										//overwrite with new record
										Reads[Amplicons[n].AmpliconID][RTI] = MakeTempRead(HeaderR1, HeaderR2, SeqR1, SeqR2, QualR1, QualR2, 
											ReadErrors, SavedBestRead.RTIErrors + RTIErrors, SavedBestRead.Frequency + 1); //increase RTI frequency

									} else {
										Reads[Amplicons[n].AmpliconID][RTI].Frequency++; //retain current record but increase frequency
										Reads[Amplicons[n].AmpliconID][RTI].RTIErrors += RTIErrors; //retain current record but increase RTIErrors
									}

								} else { //not seen before

									//bank new record
									Reads[Amplicons[n].AmpliconID][RTI] = MakeTempRead(HeaderR1, HeaderR2, SeqR1, SeqR2, QualR1, QualR2, ReadErrors, RTIErrors, 1);
								}

							} else { //?length greater than the sum of both primers
								LenDiscardedReads++;
							}

						}//?reverse primer mataches
						
						break; //if forward primer is found dont count this read again

					} //?forward primer matches

				} //end iterating over amplicons

			} //line if statement
		} //finished reading FASTQs

	} else {
		cerr << "ERROR: Unable to open FASTQ file(s)." << endl;
		return -1;
	}

	//Remove RTIs with low depth / error rate score
	RTIDepthErrorRateFilter(Reads, MinRTIDepthErrorRate);

	//Remove RTIs with an edit distance less than MinRTIEditDistance --highest frequency RTIs will be prioritised
	FilterRTIsbyEditDistance(Reads, MinRTIEditDistance);

	//print stats
	cout << "\nTotalPairedReads: " << TotalPairedReads << endl;
	cout << "N-MaskedPairedReads: " << NMaskedReads << " (" << ((float)NMaskedReads / TotalPairedReads) * 100 << "%)" << endl;
	cout << "RTIQualityDiscardedPairedReads: " << RTIQualityDiscardedReads << " (" << ((float)RTIQualityDiscardedReads / TotalPairedReads) * 100 << "%)" << endl;
	cout << "UnmatchedPrimerPairedReads: " << TotalPairedReads - (PrimerMatchedReads + RTIQualityDiscardedReads + NMaskedReads) << " (" << ((float)(TotalPairedReads - (PrimerMatchedReads + RTIQualityDiscardedReads + NMaskedReads)) / TotalPairedReads) * 100 << "%)" << endl;
	cout << "ShortInsertDiscardedPairedReads: " << LenDiscardedReads << " (" << ((float)LenDiscardedReads / TotalPairedReads) * 100 << "%)" << endl;
	cout << "Amplicon\tUsableReads\tUniqueReads\tDuplicationRate" << endl;

	//print passing records and per-amplicon stats
	for (n = 0; n < Amplicons.size(); ++n){

		for (pair<string, molecule> Read : Reads[Amplicons[n].AmpliconID]){ //RTI = molecule

			if (Read.second.PrintRead == true){

				TotalUsableMolecules++;
				AmpliconUniqueReads[Amplicons[n].AmpliconID]++;

				if (AmpliconStrand[Amplicons[n].AmpliconID] == 0){
					R1Dedupped0 << Read.second.HeaderR1 << "\012";
					R1Dedupped0 << Read.second.SeqR1 << "\012+\012";
					R1Dedupped0 << Read.second.QualR1 << "\012";

					R2Dedupped0 << Read.second.HeaderR2 << "\012";
					R2Dedupped0 << Read.second.SeqR2 << "\012+\012";
					R2Dedupped0 << Read.second.QualR2 << "\012";
				} else {
					R1Dedupped1 << Read.second.HeaderR1 << "\012";
					R1Dedupped1 << Read.second.SeqR1 << "\012+\012";
					R1Dedupped1 << Read.second.QualR1 << "\012";

					R2Dedupped1 << Read.second.HeaderR2 << "\012";
					R2Dedupped1 << Read.second.SeqR2 << "\012+\012";
					R2Dedupped1 << Read.second.QualR2 << "\012";
				}

				//AmpliconID, RTI, RTI_Frequency, RTI_ReadErrors
				StatsOut << SampleID << "\t" << Amplicons[n].AmpliconID << "\t" << AmpliconStrand[Amplicons[n].AmpliconID] << "\t" << Read.first << "\t" << Read.second.Frequency << "\t" << Read.second.ReadErrors << "\n";
			}

		}

		//print per amplicon stats
		if (AmpliconUsableReads.count(Amplicons[n].AmpliconID) == 1){ //reads associated with this amplicon
			cout << Amplicons[n].AmpliconID << "\t" << AmpliconUsableReads[Amplicons[n].AmpliconID] << "\t" << AmpliconUniqueReads[Amplicons[n].AmpliconID] << "\t" << (1 - ((float)AmpliconUniqueReads[Amplicons[n].AmpliconID] / AmpliconUsableReads[Amplicons[n].AmpliconID])) * 100 << "%" << endl;
		} else {
			cout << Amplicons[n].AmpliconID << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;
		}

	} //finish iterating over amplicons

	cout << "UniqueMolecules: " << TotalUsableMolecules << " (" << ((float)TotalUsableMolecules / TotalUsableReads) * 100 << "%)" << endl;
	cout << "DuplicationRate: " << (1 - ((float)TotalUsableMolecules / TotalUsableReads)) * 100 << "%" << endl << endl;

	//print unfiltered downsampled reads
	for (auto & Amplicon : UnfilteredReads){
		
		//randomly shuffle unfiltered paired reads
		random_shuffle(Amplicon.second.begin(), Amplicon.second.end(), RandomGenerator);

		if (AmpliconStrand[Amplicon.first] == 0){
			
			//select the first n reads giving the same depth per amplicon as filtered
			for (n = 0; n < AmpliconUniqueReads[Amplicon.first]; ++n){

				R1Trimmed0 << Amplicon.second[n].HeaderR1 << "\012";
				R1Trimmed0 << Amplicon.second[n].SeqR1 << "\012+\012";
				R1Trimmed0 << Amplicon.second[n].QualR1 << "\012";
				
				R2Trimmed0 << Amplicon.second[n].HeaderR2 << "\012";
				R2Trimmed0 << Amplicon.second[n].SeqR2 << "\012+\012";
				R2Trimmed0 << Amplicon.second[n].QualR2 << "\012";

			}

		} else {
			
			//select the first n reads giving the same depth per amplicon as filtered
			for (n = 0; n < AmpliconUniqueReads[Amplicon.first]; ++n){

				R1Trimmed1 << Amplicon.second[n].HeaderR1 << "\012";
				R1Trimmed1 << Amplicon.second[n].SeqR1 << "\012+\012";
				R1Trimmed1 << Amplicon.second[n].QualR1 << "\012";

				R2Trimmed1 << Amplicon.second[n].HeaderR2 << "\012";
				R2Trimmed1 << Amplicon.second[n].SeqR2 << "\012+\012";
				R2Trimmed1 << Amplicon.second[n].QualR2 << "\012";

			}

		}

	}

	return 0;
}