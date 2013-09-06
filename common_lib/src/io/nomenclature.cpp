#include <stdlib.h>
#include <string>
#include "nomenclature.h"
#include "generic_exception.h"
#include "tokenizer.h"

Nomenclature::Nomenclature(const std::string& data) {
	parse(data);
}

int Nomenclature::getFR1begin() const {
	return regions[1];
}

int Nomenclature::getCDR1begin() const {
	return regions[3];
}

int Nomenclature::getFR2begin() const {
	return regions[5];
}

int Nomenclature::getCDR2begin() const {
	return regions[7];
}

int Nomenclature::getFR3begin() const {
	return regions[9];
}

int Nomenclature::getFR3end() const {
	return regions[10];
}
