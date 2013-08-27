#ifndef TASK_CONFIGURATOR_H_
#define TASK_CONFIGURATOR_H_

#include "database.h"
#include "aho_corasick.h"
#include "settings.h"

class TaskConfigurator {
public:
	static void configureAndRun(const struct settings_t& settings);
private:
	static void getKmersAhoCorasick(const Database& data, AhoCorasick& dbAhoCorasick);
};


#endif /* TASK_CONFIGURATOR_H_ */
