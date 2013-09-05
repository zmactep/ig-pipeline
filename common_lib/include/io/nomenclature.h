#ifndef NOMENCLATURE_H_
#define NOMENCLATURE_H_

#include <string>

class Nomenclature {
public:
	Nomenclature(){};
	Nomenclature(const std::string& data);
	void parse(const std::string& data);
	std::string getRefName() const;
	int getFR1begin() const;
	int getCDR1begin() const;
	int getFR2begin() const;
	int getCDR2begin() const;
	int getFR3begin() const;
	int getFR3end() const;

private:
	std::string ref_name;
	int cdr1, cdr2, fr1, fr2, fr3_begin, fr3_end;
};


#endif /* NOMENCLATURE_H_ */
