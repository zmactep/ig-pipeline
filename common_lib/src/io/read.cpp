#include "read.h"

std::string Read::getName() const {
	return name;
}

std::string Read::getSeq() const {
	return seq;
}

std::string Read::getQual() const {
	return qual;
}

void Read::setName(const std::string & name) {
	this->name = name;
}

void Read::setSeq(const std::string & seq) {
	this->seq = seq;
}

void Read::setQual(const std::string & qual) {
	this->qual = qual;
}
