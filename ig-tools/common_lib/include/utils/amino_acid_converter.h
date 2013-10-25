#ifndef AMINO_ACID_CONVERTER_H_
#define AMINO_ACID_CONVERTER_H_

#include <map>
#include <vector>
#include <string>

class AminoAcidConverter {
public:
	static std::vector<std::string> convertFromAminoAcid(const std::string & amino_acids);
	static std::vector<std::string> convertFromNucleotide(const std::string & nucleo);

private:
    static std::map<int, std::vector<std::string> > create_acid_map();
    static std::map<std::string, char> create_nucleotid_map();

    static const std::map<int, std::vector<std::string> > acid2nucleotid;
    static const std::map<std::string, char> nucleotid2acid;
};


#endif /* AMINO_ACID_CONVERTER_H_ */
