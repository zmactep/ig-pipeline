#include <log4cxx/logger.h>
#include <omp.h>

#include "regions_finder.h"
#include "ssw_cpp.h"
#include "output_formatter.h"
#include "kmer_generator.h"
#include "aho_corasick.h"

extern log4cxx::LoggerPtr logger;

void RegionsFinder::operator()(const Read &r) {
  try {
    const std::string& name = r.getName();
    const std::string& sequence = r.getSeq();

    std::set<std::string *> setOfReferences2check;
    seq2index_t kmersFound = kmersAhoCorasick.search(sequence);

	std::set<std::string *, Compare> sequences;
    for (auto it = kmersFound.begin(); it != kmersFound.end(); ++it) {
    	sequences.clear();
    	data->get_sequences_for_kmer(*(it->first), sequences);
        setOfReferences2check.insert(sequences.begin(), sequences.end());
    }

    StripedSmithWaterman::Aligner aligner(match_score, mismatch_penalty, gap_opening_penalty, gap_extending_penalty);
    StripedSmithWaterman::Filter filter;
    aligner.SetReferenceSequence(sequence.c_str(), sequence.size());

    for (auto it = setOfReferences2check.begin(), et = setOfReferences2check.end(); it != et; ++it) {
      StripedSmithWaterman::Alignment alignment;
      const std::string& query = **it;
      aligner.Align(query.c_str(), filter, &alignment);

      std::string database_name;
      data->get_name_by_sequence(query, database_name);

#pragma omp critical
        {
    	  OutputFormatter::print_alignment(output_align, alignment, sequence, query, name, database_name);
    	  OutputFormatter::print_regions(output_regions, alignment, sequence, query, name, database_name);
        }
    }

  } catch (std::exception& e) {
	  LOG4CXX_DEBUG(logger, "error for read " << r.getSeq() << ":" << e.what());
  }
}
