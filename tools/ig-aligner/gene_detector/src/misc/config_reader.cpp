#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <omp.h>
#include "config_reader.h"
#include "generic_exception.h"

struct settings_t ConfigReader::read_settings() {
	std::ifstream settings_file( "./config/config.ini" );

	if (!settings_file.is_open()) {
		throw GenericException("Cannot open config file");
	}

	po::options_description desc("Options");
	desc.add_options()
			("data.v_data_file", po::value<std::string>(), "V data file")
			("data.d_data_file", po::value<std::string>(), "D data file")
			("data.j_data_file", po::value<std::string>(), "J data file")
			("data.output_file", po::value<std::string>(), "output file");
	po::variables_map vm;
	po::store(po::parse_config_file(settings_file, desc, true), vm);
	settings_file.close();
	po::notify(vm);

	struct settings_t settings;
	settings.v_data_file = vm["data.v_data_file"].as<std::string>();
	settings.d_data_file = vm["data.d_data_file"].as<std::string>();
	settings.j_data_file = vm["data.j_data_file"].as<std::string>();
	settings.output_file = vm["data.output_file"].as<std::string>();

	return settings;
}
