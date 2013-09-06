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
	std::cout << "Usage for generating train data: svm_data_generator train fasta_file nomenclature_file sliding_window_size" << std::endl;
	std::cout << "Usage for generating predict data: svm_data_generator predict fasta_file sliding_window_size" << std::endl;
}

void split_read_by_nomenclature(const std::string& r, const Nomenclature& n, std::map<int, std::string>& output) {
	//no length checks here are explicit - we better through an exception
	int offset = 0;

	for (int i = 0; i < n.getNumRegions(); ++i) {
		std::string region = r.substr(offset + n.getRegionBegin(i) - 1, n.getRegionEnd(i) - offset);
		output.insert(std::make_pair(i, region));
		offset += region.length();
	}
}

void print_svm_entry(int label, const std::string& line, const std::string& comment) {
	std::cout << label << " ";
	for (int i = 0; i < line.size(); ++i) {
		std::cout << i << ":" << (int)line[i] << " ";
	}
	std::cout << "# " << comment << std::endl;
}

void get_seq2nomenclature(const char* filename, std::map<std::string, Nomenclature>& seq2nomenclature) {
	KabatReader kabat_reader(filename);
	Nomenclature n;
	while (!kabat_reader.eof()) {
		try {
			kabat_reader >> n;
			seq2nomenclature.insert(std::make_pair(n.getRefName(), n));
		} catch (std::exception& e) {
			std::clog << "Error parsing line in kabat_reader: " << e.what() << std::endl;
		}
	}
}

//for train data
void process_read(const Read& r, const std::map<std::string, Nomenclature>& seq2nomenclature, int window_size) {
	std::map<int, std::string> regions;
	try {
		split_read_by_nomenclature(r.getSeq(), seq2nomenclature.find(r.getName())->second, regions);
	} catch (std::exception& e){
		std::clog << "Error splitting read by nomenclature: " << r.getName() << std::endl;
	}

	for (std::map<int, std::string>::const_iterator it1 = regions.begin(); it1 != regions.end(); ++it1) {
		const int label = it1->first;

		if (r.getSeq().length() > window_size) {
			std::vector<std::string> kmers = KmerGenerator::getKmers(it1->second, window_size);
			for (std::vector<std::string>::const_iterator it2 = kmers.begin(); it2 != kmers.end(); ++it2) {
				print_svm_entry(label, *it2, r.getName());
			}
		}
	}
}

//for predict data
void process_read(const Read& r, int window_size) {
	if (r.getSeq().length() > window_size) {
		std::vector<std::string> kmers = KmerGenerator::getKmers(r.getSeq(), window_size);
		for (std::vector<std::string>::const_iterator it = kmers.begin(); it != kmers.end(); ++it) {
			print_svm_entry(0, *it, r.getName());
		}
	}
}

void generate_train_date(char ** argv) {
	const int window_size = atoi(argv[4]);

	try {
		std::map<std::string, Nomenclature> seq2nomenclature;
		get_seq2nomenclature(argv[3], seq2nomenclature);

		FastaReader fasta_reader(argv[2]);
		Read r;
		while (!fasta_reader.eof()) {
			fasta_reader >> r;
			process_read(r, seq2nomenclature, window_size);
		}

	} catch (std::exception& e){
		std::cout << "Error: " << e.what() << std::endl;
	}
}

void generate_predict_date(char ** argv) {
	const int window_size = atoi(argv[3]);
	try {
		FastaReader fasta_reader(argv[2]);
		Read r;
		while (!fasta_reader.eof()) {
			fasta_reader >> r;
			process_read(r, window_size);
		}
	} catch (std::exception& e){
		std::cout << "Error: " << e.what() << std::endl;
	}
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
