#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <string>

struct settings_t {
	int max_threads, kmer_size, match_score, mismatch_penalty, gap_opening_penalty, gap_extending_penalty;
	std::string reference_file, input_file, output_align, output_regions;
};

#endif /* SETTINGS_H_ */
