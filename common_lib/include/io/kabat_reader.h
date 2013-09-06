#ifndef KABAT_READER_H_
#define KABAT_READER_H_

#include <iostream>
#include <fstream>
#include <string>
#include "generic_nomenclature.h"

class KabatReader{
public:
	KabatReader(const std::string & filename);
	virtual ~KabatReader();
	KabatReader& operator>>(GenericNomenclature &n);
	bool eof() const;
private:
	std::ifstream input;
};

#endif /* KABAT_READER_H_ */
