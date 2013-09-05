#include <iostream>
#include <exception>
#include <stdlib.h>
#include <map>
#include <string>
#include <vector>
#include "kmer_generator.h"
#include "fasta_reader.h"
#include "kabat_reader.h"
#include "nomenclature.h"
#include "read.h"

void usage(){
	std::cout << "Usage: svm_data_generator fasta_file nomenclature_file sliding_window_size" << std::endl;
}

void splitReadByNomenclature(const std::string& r, const Nomenclature& n, std::map<std::string, std::string>& output) {
	//no length checks here are explicit - we better through an exception
	int offset = 0;
	output.insert(std::make_pair("FR1", r.substr(0, n.getCDR1begin() - offset + 1)));
	offset += output["FR1"].length();
	output.insert(std::make_pair("CDR1", r.substr(n.getCDR1begin(), n.getFR2begin() - offset + 1)));
	offset += output["CDR1"].length();
	output.insert(std::make_pair("FR2", r.substr(n.getFR2begin(), n.getCDR2begin() - offset + 1)));
	offset += output["FR2"].length();
	output.insert(std::make_pair("CDR2", r.substr(n.getCDR2begin(), n.getFR3begin() - offset + 1)));
	offset += output["CDR2"].length();
	output.insert(std::make_pair("FR3", r.substr(n.getFR3begin(), n.getFR3end() - offset + 1)));
}

int main(int argc, char ** argv) {
	if (4 != argc) {
		usage();
		return 0;
	}

	const int window_size = atoi(argv[3]);

	try {
		FastaReader fasta_reader(argv[1]);
		KabatReader kabat_reader(argv[2]);

		std::map<std::string, Nomenclature> seq2nomenclature;
		Nomenclature n;
		while (!kabat_reader.eof()) {
			kabat_reader >> n;
			seq2nomenclature.insert(std::make_pair(n.getRefName(), n));
		}

		Read r;
		while (!fasta_reader.eof()) {
			fasta_reader >> r;
			std::map<std::string, std::string> regions;
			try {
				splitReadByNomenclature(r.getSeq(), n, regions);
			} catch (std::exception& e){
				std::clog << "Error processing read: " << r.getName() << std::endl;
			}

			for (std::map<std::string, std::string>::const_iterator it1 = regions.begin(); it1 != regions.end(); ++it1) {
				const std::string& label = it1->first;

				std::vector<std::string> kmers = KmerGenerator::getKmers(it1->second, window_size);

				for (std::vector<std::string>::const_iterator it2 = kmers.begin(); it2 != kmers.end(); ++it2) {
					//TODO output in libsvm format
					std::cout << label << ":" << *it2 << std::endl;
				}
			}
		}

	} catch (std::exception& e){
		std::cout << "Error: " << e.what() << std::endl;
	}
	return 0;
}
