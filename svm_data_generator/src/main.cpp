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
	std::cout << "Usage for generating train data: svm_data_generator train fasta_file GenericNomenclature_file sliding_window_size" << std::endl;
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
void process_read(const Read& r, const std::map<std::string, GenericNomenclature>& seq2GenericNomenclature, int window_size) {
	std::map<int, std::string> regions;
	std::map<std::string, GenericNomenclature>::const_iterator it;
	if (seq2GenericNomenclature.end() == (it = seq2GenericNomenclature.find(r.getName()))) {
		std::clog << "Missing nomenclature for read " << r.getName() << std::endl;
		return;
	}
	split_read_by_GenericNomenclature(r, it->second, regions);

	for (std::map<int, std::string>::const_iterator it1 = regions.begin(); it1 != regions.end(); ++it1) {
		const int label = it1->first;

		if (r.getSeq().length() > window_size) {
			std::vector<std::string> kmers = KmerGenerator::getKmers(it1->second, window_size);
			for (std::vector<std::string>::const_iterator it2 = kmers.begin(); it2 != kmers.end(); ++it2) {
				const std::string& line = *it2;
				const std::string& comment = r.getName();
				std::cout << label << " " << line << std::endl;
			}
		}
	}
}

//for predict data
void process_read(const Read& r, int window_size, std::ofstream& comments_output) {
	if (r.getSeq().length() > window_size) {
		std::vector<std::string> kmers = KmerGenerator::getKmers(r.getSeq(), window_size);
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

	try {
		std::map<std::string, GenericNomenclature> seq2nomenclature;
		get_seq2GenericNomenclature(argv[3], seq2nomenclature);

		FastaReader fasta_reader(argv[2]);
		Read r;
		while (!fasta_reader.eof()) {
			fasta_reader >> r;
			try {
				process_read(r, seq2nomenclature, window_size);
			} catch (std::exception& e) {
				std::clog << "Error processing read: " << e.what() << std::endl;
			}
		}
	} catch (std::exception& e){
		std::cout << "Error: " << e.what() << std::endl;
	}
}

void generate_predict_date(char ** argv) {
	std::ofstream output("read_names.txt");
	const int window_size = atoi(argv[3]);
	try {
		FastaReader fasta_reader(argv[2]);
		Read r;
		while (!fasta_reader.eof()) {
			fasta_reader >> r;
			try {
				process_read(r, window_size, output);
			} catch (std::exception& e) {
				std::clog << "Error processing read: " << e.what() << std::endl;
			}
		}
	} catch (std::exception& e){
		std::cout << "Error: " << e.what() << std::endl;
	}
	output.close();
}

int main(int argc, char ** argv) {
	if (5 != argc && 4 != argc) {
		usage();
		return 0;
	}

	if (!strcmp(argv[1], "train") && 5 == argc) {
		generate_train_date(argv);
	} else if (!strcmp(argv[1], "predict") && 4 == argc) {
		generate_predict_date(argv);
	} else {
		std::cout << "unknown mode: not train and not predict." << std::endl;
		usage();
	}
	return 0;
}
