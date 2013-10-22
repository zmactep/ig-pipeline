#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include "settings.h"
#include "config_reader.h"
#include "alignment.h"

void getMapForFile(std::ifstream& input, std::map<std::string, std::vector<Alignment> >& reads, std::set<std::string>& uniq_reads) {
	if (!input.is_open()) {
		return;
	}
	std::map<std::string, std::vector<Alignment> >::iterator it;
	std::string line;
	while(std::getline(input, line)) {
		try {
			Alignment a(line);
			uniq_reads.insert(a.getReadName());
			if (reads.end() != (it = reads.find(a.getReadName()))) {
				it->second.push_back(a);
			} else {
				std::vector<Alignment> alignments;
				alignments.push_back(a);
				reads.insert(std::make_pair(a.getReadName(), alignments));
			}

		} catch (std::exception & e) {
			std::cout << "Alignment error: " << e.what() << std::endl;
		}
	}

	for (it = reads.begin(); it != reads.end(); ++it) {
		std::sort(it->second.begin(), it->second.end());
		std::reverse(it->second.begin(), it->second.end());
	}
}

void fillGeneResults(const std::map<std::string, std::vector<Alignment> >& reads, const std::string& read_name, const std::string gene, std::ofstream& output) {
	auto iter = reads.find(read_name);
	output << "\t\t\"" << gene << "\": [";

	if (reads.end() != iter) {
		for (size_t i = 0; i < iter->second.size(); ++i) {
			const Alignment& a = iter->second[i];
			output << "{\"ref_name\": \"" << a.getRefName() << "\"" <<
					", \"score\": " << a.getScore() <<
					", \"mismatches\": " << a.getMismatches() <<
					", \"ref_begin\": " << a.getRefBegin() <<
					", \"ref_end\": " << a.getRefEnd() << "}";
			if (i != iter->second.size() - 1) {
				output << ", ";
			}
		}
	}
	output << "]";
}

int main() {

	std::cout << "Lets start" << std::endl;

	struct settings_t settings;

	try {
		settings = ConfigReader::read_settings();
	} catch (std::exception& e) {
		std::cout << "Config file is incorrect: " << e.what() <<std::endl;
		return 0;
	}

	std::ifstream v_data(settings.v_data_file);
	std::ifstream d_data(settings.d_data_file);
	std::ifstream j_data(settings.j_data_file);
	std::ofstream output(settings.output_file);

	std::set<std::string> uniq_reads;
	std::map<std::string, std::vector<Alignment> > v_reads, d_reads, j_reads;
	std::cout << "Reading V data..." << std::endl;
	getMapForFile(v_data, v_reads, uniq_reads);
	std::cout << "Done.\nReading D data..." << std::endl;
	getMapForFile(d_data, d_reads, uniq_reads);
	std::cout << "Done.\nReading J data..." << std::endl;
	getMapForFile(j_data, j_reads, uniq_reads);
	std::cout << "Done." << std::endl;

	v_data.close();
	d_data.close();
	j_data.close();

	std::cout << "Start processing..." << std::endl;
	int counter = 0;
	output << "{" << std::endl;
	for (auto it = uniq_reads.begin(); it != uniq_reads.end(); ++it) {
		output << "\t\"" << *it << "\": {" << std::endl;
		fillGeneResults(v_reads, *it, "V", output);
		output << "," << std::endl;
		fillGeneResults(d_reads, *it, "D", output);
		output << "," << std::endl;
		fillGeneResults(j_reads, *it, "J", output);
		output << std::endl << "\t}" << std::endl;

		if (++counter % 100) {
			std::cout << counter << " reads of " << uniq_reads.size() << " processed\r";
		}
	}
	output << std::endl << "}" << std::endl;
	output.close();

	std::cout << std::endl << "Done" << std::endl;

	return 0;
}
