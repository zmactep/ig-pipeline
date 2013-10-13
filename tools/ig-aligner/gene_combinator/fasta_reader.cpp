#include "fasta_reader.h"
#include "generic_exception.h"

bool FastaReader::open(std::string filename) {
	fp = gzopen(filename.c_str(), "r");
	if (!fp) {
		return false;
	}
	is_opened_ = true;
	seq = kseq_init(fp);
	is_eof_ = false;
	readNext();
	return true;
}

void FastaReader::readNext() {
	if(!is_open()) {
		throw GenericException("File is not opened");
	}

	if(eof()) {
		throw GenericException("End of file reached");
	}

	if (kseq_read(seq) < 0) {
		is_eof_ = true;
	}
}

FastaReader::~FastaReader(){
	close();
}

FastaReader& FastaReader::operator>>(Read &r){
	if(!is_open()) {
		throw GenericException("File is not opened");
	}

	if(eof()) {
		throw GenericException("End of file reached");
	}

	r.setName(seq->name.s);
	r.setSeq(seq->seq.s);
	if (seq->qual.s) {
		r.setQual(seq->qual.s);
	}

	readNext();

	return *this;
}

void FastaReader::reset() {
	close();
	open(filename);
}

bool FastaReader::is_open() const {
	return is_opened_;
}

bool FastaReader::eof() const{
	return is_eof_;
}

void FastaReader::close() {
	if (is_open()) {
		kseq_destroy(seq);
		gzclose(fp);
		is_opened_ = false;
		is_eof_ = false;
	}
}

