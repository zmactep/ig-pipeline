#include <string>
#include <iostream>
#include <boost/program_options.hpp>
#include <log4cxx/logger.h>
#include <log4cxx/propertyconfigurator.h>
#include <log4cxx/helpers/exception.h>
#include "config_reader.h"
#include "task_configurator.h"
#include "settings.h"

int main() {
	log4cxx::PropertyConfigurator::configure("./config/log4cxx.properties");
	log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("main"));

    LOG4CXX_INFO(logger, "Lets start");

	struct settings_t settings;

	try {
		settings = ConfigReader::read_settings();
	} catch (std::exception& e) {
		LOG4CXX_ERROR(logger, "Config file is incorrect: " << e.what());
		return 0;
	}

	clock_t start = clock();
	try {
		TaskConfigurator::configureAndRun(settings);
	} catch (std::exception& e) {
		LOG4CXX_ERROR(logger, "Config file is incorrect: " << e.what());
		return 0;
	}
	clock_t ends = clock();
	LOG4CXX_INFO(logger, "Done. Processor time spent is " << (double) (ends - start) / CLOCKS_PER_SEC << " seconds. Bye!");

    return 0;
}
