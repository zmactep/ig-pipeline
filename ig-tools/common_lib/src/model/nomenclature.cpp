#include <stdlib.h>
#include "nomenclature.h"
#include "generic_exception.h"
#include "tokenizer.h"

#define MIN_COL_NUM 11

Nomenclature::Nomenclature(const std::string& data) {
	std::vector<std::string> tokens = Tokenizer::tokenize(data);
	if (tokens.size() < MIN_COL_NUM) {
		throw GenericException("Not enough data: " + data);
	}

	this->ref_name = tokens[0];
	this->fr1 = atoi(tokens[1].c_str());
	this->cdr1 = atoi(tokens[3].c_str());
	this->fr2 = atoi(tokens[5].c_str());
	this->cdr2 = atoi(tokens[7].c_str());
	this->fr3_begin = atoi(tokens[9].c_str());
	this->fr3_end = atoi(tokens[10].c_str());
}

std::string Nomenclature::getRefName() const {
	return ref_name;
}

int Nomenclature::getFR1begin() const {
	return fr1;
}

int Nomenclature::getCDR1begin() const {
	return cdr1;
}

int Nomenclature::getFR2begin() const {
	return fr2;
}

int Nomenclature::getCDR2begin() const {
	return cdr2;
}

int Nomenclature::getFR3begin() const {
	return fr3_begin;
}

int Nomenclature::getFR3end() const {
	return fr3_end;
}
