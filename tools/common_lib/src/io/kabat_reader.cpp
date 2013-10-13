#include "kabat_reader.h"
#include "generic_exception.h"

KabatReader::KabatReader(const std::string & filename) {
	input.open(filename.c_str(), std::ifstream::in);
	if (!input.is_open()){
		throw GenericException("Cannot open file");
	}
}

KabatReader::~KabatReader() {
	input.close();
}

KabatReader& KabatReader::operator>>(GenericNomenclature &n) {
	if(!input.is_open()) {
		throw GenericException("File is not opened");
	}

	if(input.eof()) {
		throw GenericException("End of file reached");
	}

	std::string line;
	std::getline(input, line);
	n.parse(line);

	return *this;
}

bool KabatReader::eof() const {
	return input.eof();
}
