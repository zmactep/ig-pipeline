#include <vector>
#include "output_formatter.h"

void OutputFormatter::print_n_times(std::ostream& output, char c, int n) {
	for(int i = 0; i < n; ++i) {
		output << c;
	}
}

int OutputFormatter::restoreFromCigar(const std::string& ref, const std::string& query,
		std::string& out_ref, std::string& out_query, const StripedSmithWaterman::Alignment& a) {
	std::vector<char> aligned_ref, aligned_query;
	int ref_pos = 0, query_pos = 0, matches = 0;

	for (std::vector<uint32_t>::const_iterator it = a.cigar.begin(); it != a.cigar.end(); ++it) {
		int num = (*it & 0xFFFFFFF0) >> 4;
		int op_code = *it & 0x0000000F;
		switch (op_code) {
			case 0: {//match
				matches += num;
				for (int i = 0; i < num; ++i) {
					aligned_ref.push_back(ref[a.ref_begin + ref_pos++]);
					aligned_query.push_back(query[a.query_begin + query_pos++]);
				}
				break;
			}
			case 1: {//insert
				for (int i = 0; i < num; ++i) {
					aligned_ref.push_back('-');
					aligned_query.push_back(query[a.query_begin + query_pos++]);
				}
				break;
			}
			case 2: {//del
				for (int i = 0; i < num; ++i) {
					aligned_ref.push_back(ref[a.ref_begin + ref_pos++]);
					aligned_query.push_back('-');
				}
				break;
			}
			default:
				break;
		}
	}

	out_ref = std::string(aligned_ref.begin(), aligned_ref.end());
	out_query = std::string(aligned_query.begin(), aligned_query.end());
	return matches;
}

//output in blast-like style
void OutputFormatter::print_alignment(std::ostream& output, const StripedSmithWaterman::Alignment& a,
		const std::string& ref, const std::string& query,
		const std::string& ref_name, const std::string& query_name) {
	std::string aligned_query, aligned_ref;
	const int matches = restoreFromCigar(ref, query, aligned_ref, aligned_query, a);

	output << "Alignment score - " << a.sw_score << ", matches - " << matches << ", mismatches - "
			<< a.mismatches << ", cigar - " << a.cigar_string <<": " << std::endl
			<< ref_name << " (first line)  " << std::endl
			<< query_name <<" (last line) " << std::endl;


	// case when pattern's start pos is less than text one
	int text_offset = a.ref_begin - a.query_begin < 0 ? a.query_begin - a.ref_begin : 0;

	// ref string
	print_n_times(output, ' ', text_offset);
	output << ref << std::endl;
	print_n_times(output, ' ', text_offset + a.ref_begin);
	output << aligned_ref << std::endl;

	// vertical dashes
	print_n_times(output, ' ', text_offset + a.ref_begin);
	for (int i = 0; i < (int)std::min(aligned_query.length(), aligned_ref.length()); ++i) {
		aligned_query.at(i) == aligned_ref.at(i) ? output << "|" : output << "*";
	}
	output << std::endl;

	// query string
	print_n_times(output, ' ', text_offset + a.ref_begin);
	output << aligned_query << std::endl;
	print_n_times(output, ' ', a.ref_begin - a.query_begin);
	output << query << std::endl;
	output << std::endl;
}

void OutputFormatter::print_match(std::ostream& output, std::map<std::string*, std::vector<int>, Compare>& res, const std::string& name, const std::string& seq, const Database * data) {
	for (std::map<std::string*, std::vector<int>, Compare>::const_iterator it = res.begin(); it != res.end(); ++it) {
		for (std::vector<int>::const_iterator it_pos = it->second.begin(); it_pos != it->second.end(); ++it_pos) {
			std::string database_name;
			data->get_name_by_sequence(*(it->first), database_name);

			output << "Match: input sequence (first line) " << name << " matches " << std::endl
					<< "sequence from database (2nd line) " << database_name << std::endl;

			output << seq << std::endl;
			print_n_times(output, ' ', *it_pos);
			print_n_times(output, '|', it->first->length());
			output << std::endl;
			print_n_times(output, ' ', *it_pos);
			output << *(it->first) << std::endl;
			output << std::endl;
		}
	}
}

void OutputFormatter::print_regions(std::ostream& output, const StripedSmithWaterman::Alignment& a,
		const std::string& ref, const std::string& query,
		const std::string& ref_name, const std::string& query_name) {
	output << ref_name << "\t" << query_name << "\t" << a.ref_begin << "\t" << a.ref_end << "\t" << a.query_begin << "\t" << a.query_end << "\t" << a.mismatches << "\t"<< a.sw_score << std::endl;
}
