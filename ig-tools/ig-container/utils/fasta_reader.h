#ifndef FASTA_READER_H_
#define FASTA_READER_H_

#include <zlib.h>
#include <string>
#include "kseq.h"
#include "read.h"

KSEQ_INIT(gzFile, gzread)

class FastaReader{
public:
	FastaReader(const std::string & filename) : filename(filename) {
		open(filename);
	}

	virtual ~FastaReader();
	FastaReader& operator>>(Read &r);
	void reset();
	bool is_opened() const;
	bool is_eof() const;
private:
	bool open(std::string filename);
	void readNext();
	void close();
	std::string filename;
	bool is_opened_;
	bool is_eof_;
	gzFile fp;
	kseq_t* seq;
};


#endif /* FASTA_READER_H_ */
