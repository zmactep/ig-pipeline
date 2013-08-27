#include <stdlib.h>
#include <vector>
#include "generic_exception.h"
#include "alignment.h"
#include "tokenizer.h"

#define TOKENS_NUM 8

Alignment::Alignment(const std::string& data) {
	std::vector<std::string> tokens = Tokenizer::tokenize(data);
	if (TOKENS_NUM != tokens.size()) {
		throw GenericException("Not enough data: " + data);
	}

	this->read_name = tokens[0];
	this->ref_name = tokens[1];
	this->read_begin = atoi(tokens[2].c_str());
	this->read_end = atoi(tokens[3].c_str());
	this->ref_begin = atoi(tokens[4].c_str());
	this->ref_end = atoi(tokens[5].c_str());
	this->mismatches = atoi(tokens[6].c_str());
	this->score = atoi(tokens[7].c_str());
}

bool Alignment::operator<(const Alignment& rhs) const {
	return this->score < rhs.score;
}

bool Alignment::operator>(const Alignment& rhs) const {
	return this->score < rhs.score;
}

std::string Alignment::getRefName() const {
	return ref_name;
}

std::string Alignment::getReadName() const {
	return read_name;
}

int Alignment::getReadBegin() const {
	return read_begin;
}

int Alignment::getReadEnd() const {
	return read_end;
}

int Alignment::getRefBegin() const {
	return ref_begin;
}

int Alignment::getRefEnd() const {
	return ref_end;
}

int Alignment::getMismatches() const {
	return mismatches;
}

int Alignment::getScore() const {
	return score;
}
