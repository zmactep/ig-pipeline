#ifndef GENERIC_NOMENCLATURE_H_
#define GENERIC_NOMENCLATURE_H_

#include <string>
#include <vector>

class GenericNomenclature {
public:
	GenericNomenclature(){};
	GenericNomenclature(const std::string& data);
	~GenericNomenclature(){};
	void parse(const std::string& data);

	std::string getRefName() const;
	int getNumRegions() const;
	int getRegionBegin(int number) const;
	int getRegionEnd(int number) const;
protected:
	std::string ref_name;
	std::vector<int> regions;
};


#endif /* GENERIC_NOMENCLATURE_H_ */
