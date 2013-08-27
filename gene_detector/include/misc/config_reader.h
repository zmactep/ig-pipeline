#ifndef CONFIG_READER_H_
#define CONFIG_READER_H_

#include <boost/program_options.hpp>
#include "settings.h"

namespace po = boost::program_options;

class ConfigReader {
public:
	static struct settings_t read_settings();
};

#endif /* CONFIG_READER_H_ */
