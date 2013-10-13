#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include "settings.h"
#include "config_reader.h"
#include "alignment.h"

void getMapForFile(std::ifstream& input, std::map<std::string, Alignment>& reads, std::set<std::string>& uniq_reads) {
	if (!input.is_open()) {
		return;
	}
	std::map<std::string, Alignment>::iterator it;
	std::string line;
	while(std::getline(input, line)) {
		try {
			Alignment a(line);
			uniq_reads.insert(a.getReadName());
			if (reads.end() != (it = reads.find(a.getReadName()))) {
				if (it->second < a) {
					reads.erase(it);
					reads.insert(std::make_pair(a.getReadName(), a)); //save alignment with max score
				}
			} else {
				reads.insert(std::make_pair(a.getReadName(), a));
			}

		} catch (std::exception & e) {
			std::cout << "Alignment error: " << e.what() << std::endl;
		}
	}
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
	std::map<std::string, Alignment> v_reads, d_reads, j_reads;
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
	for (auto it = uniq_reads.begin(); it != uniq_reads.end(); ++it) {
		output << *it << "\t";
		auto iter = v_reads.find(*it);
		if (v_reads.end() == iter) {
			output << "-\t-\t";
		} else {
			output << iter->second.getReadBegin() << "\t" << iter->second.getReadEnd() << "\t";
		}

		iter = d_reads.find(*it);
		if (d_reads.end() == iter) {
			output << "-\t-\t";
		} else {
			output << iter->second.getReadBegin() << "\t" << iter->second.getReadEnd() << "\t";
		}

		iter = j_reads.find(*it);
		if (j_reads.end() == iter) {
			output << "-\t-\t";
		} else {
			output << iter->second.getReadBegin() << "\t" << iter->second.getReadEnd() << "\t";
		}
		output << std::endl;
		if (++counter % 100) {
			std::cout << counter << " reads of " << uniq_reads.size() << " processed\r";
		}
	}

	output.close();
	std::cout << std::endl << "Done" << std::endl;

	return 0;
}
