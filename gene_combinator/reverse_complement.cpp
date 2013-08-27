#include <map>
#include <vector>
#include <algorithm>
#include "reverse_complement.h"
#include "generic_exception.h"

std::string ReverseComplement::getRevc(const std::string& data) {
	std::map<char, char> reverse;
	reverse['C'] = 'G';
	reverse['G'] = 'C';
	reverse['T'] = 'A';
	reverse['A'] = 'T';
	reverse['N'] = 'N';

	std::vector<char> res;
	for(int i = 0; i < (int) data.length(); ++i) {
		std::map<char, char>::const_iterator it = reverse.find(data[i]);
		if (reverse.end() == it) {
			throw GenericException("Unknown nucleotid");
		}
		res.push_back(it->second);
	}

	std::reverse(res.begin(), res.end());
	return std::string(res.begin(), res.end());
}
