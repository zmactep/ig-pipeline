#ifndef OUTPUT_FORMATTER_H_
#define OUTPUT_FORMATTER_H_

#include <string>
#include <map>
#include <iostream>
#include "ssw_cpp.h"
#include "comparator.h"
#include "database.h"

class OutputFormatter {
public:
	static void print_alignment(std::ostream& output, const StripedSmithWaterman::Alignment& a,
			const std::string& ref, const std::string& query, const std::string& ref_name,
			const std::string& query_name);

	static void print_regions(std::ostream& output, const StripedSmithWaterman::Alignment& a,
			const std::string& ref, const std::string& query, const std::string& ref_name,
			const std::string& query_name);

	static void print_match(std::ostream& output, std::map<std::string*, std::vector<int>, Compare>& res,
			const std::string& name, const std::string& seq, const Database * data);

private:
	static void print_n_times(std::ostream& output, char c, int n);
	static int restoreFromCigar(const std::string& ref, const std::string& query, std::string& out_ref,
			std::string& out_query, const StripedSmithWaterman::Alignment& a);
};


#endif /* OUTPUT_FORMATTER_H_ */
