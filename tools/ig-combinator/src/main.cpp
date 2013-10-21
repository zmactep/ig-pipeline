#include <iostream>
#include <fstream>
#include <sstream>
#include "read.h"
#include "fasta_reader.h"
#include "reverse_complement.h"

#define FORWARD_READING_FRAMES_COUNT 3
#define FORWARD_AND_REVC_READING_FRAMES_COUNT 6
#define FASTA_LINE_LEN 60

void print_output(std::ostream& output, const std::string& name, const std::string& sequence) {
	output << ">" << name << std::endl;

	const int iterations = sequence.length() % FASTA_LINE_LEN == 0
			? sequence.length() / FASTA_LINE_LEN
			: sequence.length() / FASTA_LINE_LEN + 1;
	for (int i = 0; i < iterations; ++i) {
		output << sequence.substr(i * FASTA_LINE_LEN, FASTA_LINE_LEN) << std::endl;
	}

	output << std::endl;
}

void print_usage() {
	std::cout << "Usage: gene_combinator first_genes.fasta last_genes.fasta output.fasta revc{0,1}" << std::endl;
	std::cout << "revc - do we need to take reverse complements into account? Either 0 (we do not) or 1 (we do)" << std::endl;
}

int main(int argc, char** argv){
	if (5 != argc) {
		print_usage();
		return 0;
	}
	std::cout << "Welcome to gene combinator." << std::endl;

	FastaReader first_gene(argv[1]);
	FastaReader last_gene(argv[2]);
	std::ofstream output(argv[3]);
	const int reading_frames = atoi(argv[4]) == 0 ? FORWARD_READING_FRAMES_COUNT : FORWARD_AND_REVC_READING_FRAMES_COUNT;

	if (!first_gene.is_open() || !last_gene.is_open() || !output.is_open()) {
		std::cout << "Error opening files." << std::endl;
	}

	int counter = 0;
	while (!first_gene.eof()) {
		Read r1;
		first_gene >> r1;

		while (!last_gene.eof()) {
			Read r2;
			last_gene >> r2;

			for (int i = 0; i < FORWARD_READING_FRAMES_COUNT; ++i) {
				std::ostringstream s;
				s << r1.getName() << "_" << r2.getName() << "{" << i << "}";
				const std::string name(s.str());
				const std::string sequence = r1.getSeq() + r2.getSeq().substr(i);
				print_output(output, name, sequence);
			}

			for (int i = FORWARD_READING_FRAMES_COUNT; i < reading_frames; ++i) {
				std::ostringstream s;
				s << r1.getName() << "_" << r2.getName() << "{" << i << "}";
				const std::string name(s.str());
				const std::string sequence = r1.getSeq() + ReverseComplement::getRevc(r2.getSeq()).substr(i - FORWARD_READING_FRAMES_COUNT);
				print_output(output, name, sequence);
			}
		}

		last_gene.reset();
		std::cout << "processed: " << ++counter << "\r";
 	}

	output.close();
	std::cout << std::endl << "Goodbye" << std::endl;
	return 0;
}
