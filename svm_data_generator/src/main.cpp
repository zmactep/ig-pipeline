#include <iostream>
#include <exception>
#include <stdlib.h>
#include <map>
#include <string>
#include <vector>
#include "kmer_generator.h"
#include "fasta_reader.h"
#include "kabat_reader.h"
#include "generic_nomenclature.h"
#include "read.h"

void usage(){
	std::cout << "Usage for generating train data: svm_data_generator train fasta_file GenericNomenclature_file sliding_window_size no_region_borders_flag" << std::endl;
	std::cout << "Usage for generating predict data: svm_data_generator predict fasta_file sliding_window_size" << std::endl;
}

//TODO move to separate lib
void split_read_by_GenericNomenclature(const Read& r, const GenericNomenclature& n, std::map<int, std::string>& output) {
	int i = 0;
	try {
		for (; i < n.getNumRegions(); ++i) {
			std::string region = r.getSeq().substr(n.getRegionBegin(i) - 1, n.getRegionEnd(i) - n.getRegionBegin(i) + 1);
			output.insert(std::make_pair(i, region));
		}
	} catch (std::exception& e){
		const int start = n.getRegionBegin(i) - 1;
		const int size = n.getRegionEnd(i) - n.getRegionBegin(i) + 1;
			std::clog << "Error splitting read by GenericNomenclature: " << r.getName() <<
					": substring index out of range: start = " << start << "; size = " << size << std::endl;
	}
}

void print_svm_entry(int label, const std::string& line, const std::string& comment) {
}

void get_seq2GenericNomenclature(const char* filename, std::map<std::string, GenericNomenclature>& seq2GenericNomenclature) {
	KabatReader kabat_reader(filename);
	GenericNomenclature n;
	while (!kabat_reader.eof()) {
		try {
			kabat_reader >> n;
			seq2GenericNomenclature.insert(std::make_pair(n.getRefName(), n));
		} catch (std::exception& e) {
			std::clog << "Error parsing line in kabat_reader: " << e.what() << std::endl;
		}
	}
}

//for train data
void process_read(const Read& r, const std::map<std::string, GenericNomenclature>& seq2GenericNomenclature, int window_size, int no_borders) {
	std::map<int, std::string> regions;
	std::map<std::string, GenericNomenclature>::const_iterator it;
	if (seq2GenericNomenclature.end() == (it = seq2GenericNomenclature.find(r.getName()))) {
		std::clog << "Missing nomenclature for read " << r.getName() << std::endl;
		return;
	}

	const GenericNomenclature& nomenclature = it->second;
	const std::string cap(window_size / 2, 'B'); //'B' is not an amino acid nor a nucleotide
	const std::string capped_seq = cap + r.getSeq() + cap; //add caps to sequence to deal with sliding window near the edges
	if (!no_borders) {
		split_read_by_GenericNomenclature(r, nomenclature, regions);

		for (std::map<int, std::string>::const_iterator it1 = regions.begin(); it1 != regions.end(); ++it1) {
			const int label = it1->first;

			if (capped_seq.length() > window_size) {
				std::vector<std::string> kmers = KmerGenerator::getKmers(it1->second, window_size);
				for (std::vector<std::string>::const_iterator it2 = kmers.begin(); it2 != kmers.end(); ++it2) {
					const std::string& line = *it2;
					std::cout << label << " " << line << std::endl;
				}
			}
		}
	} else {
		if (capped_seq.length() > window_size) {
			std::vector<std::string> kmers = KmerGenerator::getKmers(capped_seq, window_size);
			int current_kmer = 0;
			for (int current_region = 0; current_region < nomenclature.getNumRegions(); ++current_region) {
				for (int j = nomenclature.getRegionBegin(current_region); j <= nomenclature.getRegionEnd(current_region); ++j) {
					std::cout << current_region << " " << kmers[current_kmer++] << std::endl;
				}
			}
		}
	}
}

//for predict data
void process_read(const Read& r, int window_size, std::ofstream& comments_output) {
	const std::string cap(window_size / 2, 'B'); //'B' is not an amino acid nor a nucleotide
	const std::string capped_seq = cap + r.getSeq() + cap; //add caps to sequence to deal with sliding window near the edges

	if (capped_seq.length() > window_size) {
		std::vector<std::string> kmers = KmerGenerator::getKmers(capped_seq, window_size);
		for (std::vector<std::string>::const_iterator it = kmers.begin(); it != kmers.end(); ++it) {
			const std::string& line = *it;
			const std::string& comment = r.getName();
			std::cout << 0 << " " << line << std::endl;
			comments_output << comment << std::endl;
		}
	}
}

void generate_train_date(char ** argv) {
	const int window_size = atoi(argv[4]);
	const int no_borders = atoi(argv[5]);

	try {
		std::map<std::string, GenericNomenclature> seq2nomenclature;
		get_seq2GenericNomenclature(argv[3], seq2nomenclature);

		FastaReader fasta_reader(argv[2]);
		Read r;
		while (!fasta_reader.eof()) {
			fasta_reader >> r;
			try {
				process_read(r, seq2nomenclature, window_size, no_borders);
			} catch (std::exception& e) {
				std::clog << "Error processing read: " << e.what() << std::endl;
			}
		}
	} catch (std::exception& e){
		std::cout << "Error: " << e.what() << std::endl;
	}
}

void generate_predict_date(char ** argv) {
	std::ofstream comments_output("read_names.txt");
	const int window_size = atoi(argv[3]);
	try {
		FastaReader fasta_reader(argv[2]);
		Read r;
		while (!fasta_reader.eof()) {
			fasta_reader >> r;
			try {
				process_read(r, window_size, comments_output);
			} catch (std::exception& e) {
				std::clog << "Error processing read: " << e.what() << std::endl;
			}
		}
	} catch (std::exception& e){
		std::cout << "Error: " << e.what() << std::endl;
	}
	comments_output.close();
}

int main(int argc, char ** argv) {
	if (6 != argc && 4 != argc) {
		usage();
		return 0;
	}

	if (!strcmp(argv[1], "train") && 6 == argc) {
		generate_train_date(argv);
	} else if (!strcmp(argv[1], "predict") && 4 == argc) {
		generate_predict_date(argv);
	} else {
		std::cout << "unknown mode: not train and not predict." << std::endl;
		usage();
	}
	return 0;
}
