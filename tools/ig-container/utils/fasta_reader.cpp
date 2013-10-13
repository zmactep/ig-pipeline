#include "fasta_reader.h"

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
	if(!is_opened()) {
		//throw GenericException("File is not opened");
	}

	if(is_eof()) {
		//throw GenericException("End of file reached");
	}

	if (kseq_read(seq) < 0) {
		is_eof_ = true;
	}
}

FastaReader::~FastaReader(){
	close();
}

FastaReader& FastaReader::operator>>(Read &r){
	if(!is_opened()) {
		//throw GenericException("File is not opened");
	}

	if(is_eof()) {
		//throw GenericException("End of file reached");
	}

	r.name = seq->name.s;
	r.seq = seq->seq.s;
	if (seq->qual.s) {
		r.qual = seq->qual.s;
	}

	readNext();

	return *this;
}

void FastaReader::reset() {
	close();
	open(filename);
}

bool FastaReader::is_opened() const {
	return is_opened_;
}

bool FastaReader::is_eof() const{
	return is_eof_;
}

void FastaReader::close() {
	if (is_opened()) {
		kseq_destroy(seq);
		gzclose(fp);
		is_opened_ = false;
		is_eof_ = false;
	}
}

