#include "amino_acid_converter.h"
#include "generic_exception.h"

std::vector<std::string> AminoAcidConverter::convertFromAminoAcid(const std::string & amino_acids) {
	std::vector<std::string> result;
	for (unsigned int i = 0; i < amino_acids.size(); ++i) {
		char n = amino_acids[i];

		std::map<int, std::vector<std::string> >::const_iterator it;
		if (acid2nucleotid.end() == (it = acid2nucleotid.find(n))) {
			throw GenericException("Unknown amino acid");
		}

		std::vector<std::string> tmp, suffixes(it->second.begin(), it->second.end());

		if (result.empty()) {
			result.swap(suffixes);
			continue;
		}

		for (std::vector<std::string>::const_iterator it1 = result.begin(); it1 != result.end(); ++it1) {
			for (std::vector<std::string>::const_iterator it2 = suffixes.begin(); it2 != suffixes.end(); ++it2) {
				tmp.push_back(*it1 + *it2);
			}
		}

		result.swap(tmp);
	}

	return result;
}

std::vector<std::string> AminoAcidConverter::convertFromNucleotide(const std::string & nucleo) {
	std::vector<std::string> result;
	int i = 0;

	for (int i = 0; i < 3; ++i) {
		std::vector<char> res;
		int j = i;
		while(j <= nucleo.size() - 3) {
			std::map<std::string, char>::const_iterator it = nucleotid2acid.find(nucleo.substr(j, 3));
			if (nucleotid2acid.end() == it) {
				throw GenericException("Unknown nucleotid sequence");
			}
			const char c = it->second;
			res.push_back(c);
			j += 3;
		}
		result.push_back(std::string(res.begin(), res.end()));
	}
	return result;
}

const std::map<int, std::vector<std::string> > AminoAcidConverter::acid2nucleotid = AminoAcidConverter::create_acid_map();
const std::map<std::string, char> AminoAcidConverter::nucleotid2acid = AminoAcidConverter::create_nucleotid_map();

std::map<int, std::vector<std::string> > AminoAcidConverter::create_acid_map() {
	std::map<int, std::vector<std::string> > dict;

	std::string alanine[] 		= 	{"GCT", "GCC", "GCA", "GCG"},
	arginine[] 		= 	{"AGA", "AGG", "CGT", "CGC", "CGA", "CGG"},
	asparagine[] 	= 	{"AAT", "AAC"},
	aspartic_acid[] = 	{"GAT", "GAC"},
	cysteine[]		= 	{"TGT", "TGC"},
	glutamic_acid[] = 	{"GAA", "GAG"},
	glutamine[] 	= 	{"CAA", "CAG"},
	glycine[] 		= 	{"GGT", "GGC", "GGA", "GGG"},
	histidine[] 	= 	{"CAT", "CAC"},
	isoleucine[] 	= 	{"ATT", "ATC", "ATA"},

	leucine[] 		= 	{"CTT", "CTC", "CTA", "CTG", "TTG", "TTA"},
	lysine[] 		= 	{"AAA", "AAG"},
	methionine[] 	= 	{"ATG"},
	phenylalanine[] = 	{"TTT", "TTC"},
	proline[] 		= 	{"CCT", "CCC", "CCA", "CCG"},
	serine[] 		= 	{"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"},
	threonine[] 	= 	{"ACT", "ACC", "ACA", "ACG"},
	tryptophan[] 	= 	{"TGG"},
	tyrosine[] 		= 	{"TAT", "TAC"},
	valine[] 		= 	{"GTT", "GTC", "GTA", "GTG"};

	dict.insert(std::make_pair('A', std::vector<std::string>(alanine, alanine + (sizeof(alanine) / sizeof(std::string)))));
	dict.insert(std::make_pair('R', std::vector<std::string>(arginine, arginine + (sizeof(arginine) / sizeof(std::string)))));
	dict.insert(std::make_pair('N', std::vector<std::string>(asparagine, asparagine + (sizeof(asparagine) / sizeof(std::string)))));
	dict.insert(std::make_pair('D', std::vector<std::string>(aspartic_acid, aspartic_acid + (sizeof(aspartic_acid) / sizeof(std::string)))));
	dict.insert(std::make_pair('C', std::vector<std::string>(cysteine, cysteine + (sizeof(cysteine) / sizeof(std::string)))));
	dict.insert(std::make_pair('E', std::vector<std::string>(glutamic_acid, glutamic_acid + (sizeof(glutamic_acid) / sizeof(std::string)))));
	dict.insert(std::make_pair('Q', std::vector<std::string>(glutamine, glutamine + (sizeof(glutamine) / sizeof(std::string)))));
	dict.insert(std::make_pair('G', std::vector<std::string>(glycine, glycine + (sizeof(glycine) / sizeof(std::string)))));
	dict.insert(std::make_pair('H', std::vector<std::string>(histidine, histidine + (sizeof(histidine) / sizeof(std::string)))));
	dict.insert(std::make_pair('I', std::vector<std::string>(isoleucine, isoleucine + (sizeof(isoleucine) / sizeof(std::string)))));

	dict.insert(std::make_pair('L', std::vector<std::string>(leucine, leucine + (sizeof(leucine) / sizeof(std::string)))));
	dict.insert(std::make_pair('K', std::vector<std::string>(lysine, lysine + (sizeof(lysine) / sizeof(std::string)))));
	dict.insert(std::make_pair('M', std::vector<std::string>(methionine, methionine + (sizeof(methionine) / sizeof(std::string)))));
	dict.insert(std::make_pair('F', std::vector<std::string>(phenylalanine, phenylalanine + (sizeof(phenylalanine) / sizeof(std::string)))));
	dict.insert(std::make_pair('P', std::vector<std::string>(proline, proline + (sizeof(proline) / sizeof(std::string)))));
	dict.insert(std::make_pair('S', std::vector<std::string>(serine, serine + (sizeof(serine) / sizeof(std::string)))));
	dict.insert(std::make_pair('T', std::vector<std::string>(threonine, threonine + (sizeof(threonine) / sizeof(std::string)))));
	dict.insert(std::make_pair('W', std::vector<std::string>(tryptophan, tryptophan + (sizeof(tryptophan) / sizeof(std::string)))));
	dict.insert(std::make_pair('Y', std::vector<std::string>(tyrosine, tyrosine + (sizeof(tyrosine) / sizeof(std::string)))));
	dict.insert(std::make_pair('V', std::vector<std::string>(valine, valine + (sizeof(valine) / sizeof(std::string)))));

	return dict;
}

std::map<std::string, char> AminoAcidConverter::create_nucleotid_map() {
	std::map<std::string, char> dict;
	dict["TTT"] = 'F';
	dict["TTC"] = 'F';
	dict["TTA"] = 'L';
	dict["TTG"] = 'L';
	dict["CTT"] = 'L';
	dict["CTC"] = 'L';
	dict["CTA"] = 'L';
	dict["CTG"] = 'L';
	dict["ATT"] = 'I';
	dict["ATC"] = 'I';
	dict["ATA"] = 'I';
	dict["ATG"] = 'M';
	dict["GTT"] = 'V';
	dict["GTC"] = 'V';
	dict["GTA"] = 'V';
	dict["GTG"] = 'V';
	dict["TCT"] = 'S';
	dict["TCC"] = 'S';
	dict["TCA"] = 'S';
	dict["TCG"] = 'S';
	dict["CCT"] = 'P';
	dict["CCC"] = 'P';
	dict["CCA"] = 'P';
	dict["CCG"] = 'P';
	dict["ACT"] = 'T';
	dict["ACC"] = 'T';
	dict["ACA"] = 'T';
	dict["ACG"] = 'T';
	dict["GCT"] = 'A';
	dict["GCC"] = 'A';
	dict["GCA"] = 'A';
	dict["GCG"] = 'A';
	dict["TAT"] = 'Y';
	dict["TAC"] = 'Y';
	dict["TAA"] = 'z'; //stop codon
	dict["TAG"] = 'z'; //stop codon
	dict["CAT"] = 'H';
	dict["CAC"] = 'H';
	dict["CAA"] = 'Q';
	dict["CAG"] = 'Q';
	dict["AAT"] = 'N';
	dict["AAC"] = 'N';
	dict["AAA"] = 'K';
	dict["AAG"] = 'K';
	dict["GAT"] = 'D';
	dict["GAC"] = 'D';
	dict["GAA"] = 'E';
	dict["GAG"] = 'E';
	dict["TGT"] = 'C';
	dict["TGC"] = 'C';
	dict["TGA"] = 'z'; //stop codon
	dict["TGG"] = 'W';
	dict["CGT"] = 'R';
	dict["CGC"] = 'R';
	dict["CGA"] = 'R';
	dict["CGG"] = 'R';
	dict["AGT"] = 'S';
	dict["AGC"] = 'S';
	dict["AGA"] = 'R';
	dict["AGG"] = 'R';
	dict["GGT"] = 'G';
	dict["GGC"] = 'G';
	dict["GGA"] = 'G';
	dict["GGG"] = 'G';
	return dict;
}
