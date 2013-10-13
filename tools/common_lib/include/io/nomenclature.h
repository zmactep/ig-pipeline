#ifndef NOMENCLATURE_H_
#define NOMENCLATURE_H_

#include <string>
#include "generic_nomenclature.h"

//kabat nomenclature
class Nomenclature : public GenericNomenclature {
public:
	Nomenclature(){};
	Nomenclature(const std::string& data);
	~Nomenclature(){};

	int getFR1begin() const;
	int getCDR1begin() const;
	int getFR2begin() const;
	int getCDR2begin() const;
	int getFR3begin() const;
	int getFR3end() const;
};


#endif /* NOMENCLATURE_H_ */
