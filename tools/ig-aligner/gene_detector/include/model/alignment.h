#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

class Alignment {
public:
	Alignment(const std::string& data);
	bool operator<(const Alignment& rhs) const;
	bool operator>(const Alignment& rhs) const;
	std::string getRefName() const;
	std::string getReadName() const;
	int getReadBegin() const;
	int getReadEnd() const;
	int getRefBegin() const;
	int getRefEnd() const;
	int getMismatches() const;
	int getScore() const;

private:
	std::string ref_name, read_name;
	int read_begin, read_end, ref_begin, ref_end, mismatches, score;
};


#endif /* ALIGNMENT_H_ */
