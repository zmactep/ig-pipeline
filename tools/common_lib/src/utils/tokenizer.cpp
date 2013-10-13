#include <algorithm>
#include <sstream>
#include <iostream>
#include <iterator>
#include "tokenizer.h"

std::vector<std::string> Tokenizer::tokenize(const std::string& line) {
	std::vector<std::string> tokens;
    std::istringstream iss(line);
	std::copy(std::istream_iterator<std::string>(iss),
			std::istream_iterator<std::string>(),
			std::back_inserter<std::vector<std::string> >(tokens));
	return tokens;
}

