#ifndef KMER_GENERATOR_H_
#define KMER_GENERATOR_H_

#include <vector>
#include <string>
#include "read.h"

class KmerGenerator {
public:
	static std::vector<std::string> getKmers(const std::string& read, int k);
	static std::vector<std::string> getKmers(const Read& read, int k);
};


#endif /* KMER_GENERATOR_H_ */
