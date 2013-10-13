#ifndef READ_H_
#define READ_H_

#include <string>

class Read {
public:
	Read(const std::string & name, const std::string & seq, const std::string & qual) : name(name), seq(seq), qual(qual) {};
	Read(const std::string & name, const std::string & seq) : name(name), seq(seq){};
	Read(){};

	std::string getName() const;
	std::string getSeq() const;
	std::string getQual() const;
	void setName(const std::string & name);
	void setSeq(const std::string & seq);
	void setQual(const std::string & qual);

private:
	std::string name, seq, qual;
};


#endif /* READ_H_ */
