#ifndef REGIONS_FINDER_H_
#define REGIONS_FINDER_H_

#include <iostream>
#include "database.h"
#include "settings.h"
#include "aho_corasick.h"

class RegionsFinder {
public:
	RegionsFinder(const Database * data, std::ostream& output_align, std::ostream& output_regions,
			const AhoCorasick &kmersAhoCorasick, const struct settings_t& settings)
      : data(data), output_align(output_align), output_regions(output_regions),
        kmersAhoCorasick(kmersAhoCorasick), kmer_size(settings.kmer_size),
        match_score(settings.match_score), mismatch_penalty(settings.mismatch_penalty),
        gap_opening_penalty(settings.gap_opening_penalty), gap_extending_penalty(settings.gap_extending_penalty){
		output_regions << "read1\tread2\tread1.begin\tread1.end\tread2.begin\tread2.end\tmismatches\tscore" << std::endl;
	};

	void operator()(const Read &r);

private:
  const Database * data;
  std::ostream& output_align;
  std::ostream& output_regions;
  AhoCorasick kmersAhoCorasick;
  const int kmer_size;
  const int match_score;
  const int mismatch_penalty;
  const int gap_opening_penalty;
  const int gap_extending_penalty;
};

#endif /* REGIONS_FINDER_H_ */
