#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <omp.h>
#include "config_reader.h"

struct settings_t ConfigReader::read_settings() {
	std::ifstream settings_file( "./config/config.ini" );

	po::options_description desc("Options");
	desc.add_options()
			("data.reference_file", po::value<std::string>(), "reference file with igs")
			("data.input_file", po::value<std::string>(), "input fasta files")
			("data.output_align", po::value<std::string>(), "output aligned data")
			("data.output_regions", po::value<std::string>(), "output regions data")
			("settings.max_nthreads", po::value<std::string>(), "max threads")
			("settings.kmer_size", po::value<std::string>(), "kmers size")
			("settings.match_score", po::value<std::string>(), "match score")
			("settings.mismatch_penalty", po::value<std::string>(), "mismatch penalty")
			("settings.gap_opening_penalty", po::value<std::string>(), "gap opening penalty")
			("settings.gap_extending_penalty", po::value<std::string>(), "gap extending penalty");
	po::variables_map vm;
	po::store(po::parse_config_file(settings_file, desc, true), vm);
	settings_file.close();
	po::notify(vm);

	struct settings_t settings;
	settings.reference_file = vm["data.reference_file"].as<std::string>();
	settings.input_file = vm["data.input_file"].as<std::string>();
	settings.output_align = vm["data.output_align"].as<std::string>();
	settings.output_regions = vm["data.output_regions"].as<std::string>();
	settings.max_threads = std::min(atoi((vm["settings.max_nthreads"].as<std::string>()).c_str()), omp_get_max_threads());
	settings.kmer_size = atoi((vm["settings.kmer_size"].as<std::string>()).c_str());
	settings.match_score = atoi((vm["settings.match_score"].as<std::string>()).c_str());
	settings.mismatch_penalty = atoi((vm["settings.mismatch_penalty"].as<std::string>()).c_str());
	settings.gap_opening_penalty = atoi((vm["settings.gap_opening_penalty"].as<std::string>()).c_str());
	settings.gap_extending_penalty = atoi((vm["settings.gap_extending_penalty"].as<std::string>()).c_str());

	return settings;
}
