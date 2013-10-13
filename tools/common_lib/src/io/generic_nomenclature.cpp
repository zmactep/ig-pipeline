#include <stdlib.h>
#include <iostream>
#include <string>
#include "nomenclature.h"
#include "generic_exception.h"
#include "tokenizer.h"
#define MIN_COL_NUM 11

GenericNomenclature::GenericNomenclature(const std::string& data) {
	parse(data);
}

std::string GenericNomenclature::getRefName() const {
	return ref_name;
}

void GenericNomenclature::parse(const std::string& data) {
	if (data.empty()) {
		throw GenericException("Empty line in GenericNomenclature::parse");
	}
	std::vector<std::string> tokens = Tokenizer::tokenize(data);
	this->ref_name = tokens[0];
	regions.clear();;

	try {
		for (unsigned int i = 1; i < tokens.size(); ++i) {
			regions.push_back(std::stoi(tokens[i]));
		}
	} catch (std::exception& e) {
		//std::clog << "Parse error: " << data << std::endl;
	}
}

int GenericNomenclature::getNumRegions() const {
	return regions.size() / 2;
}

int GenericNomenclature::getRegionBegin(int number) const {
	if (2 * number > this->regions.size()) {
		throw GenericException("Region begin index out of range in: " + ref_name);
	}
	return regions[2 * number];
}

int GenericNomenclature::getRegionEnd(int number) const {
	if (2 * number > this->regions.size()) {
		throw GenericException("Region end index out of range in: " + ref_name);
	}
	return regions[2 * number + 1];
}
