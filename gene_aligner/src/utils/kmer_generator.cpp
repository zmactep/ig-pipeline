#include "kmer_generator.h"

std::vector<std::string> KmerGenerator::getKmers(const std::string& read, int k) {
	std::vector<std::string> result;
	for (int i = 0; i < read.length() - k + 1; ++i) {
		result.push_back(read.substr(i, k));
	}
	return result;
}

std::vector<std::string> KmerGenerator::getKmers(const Read& read, int k) {
	return KmerGenerator::getKmers(read.getSeq(), k);
}

