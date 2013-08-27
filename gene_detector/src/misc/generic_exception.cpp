#include "generic_exception.h"

GenericException::GenericException(const std::string& m) {
	this->message = m;
};

GenericException::~GenericException() throw (){};

const char* GenericException::what() const throw () {
	return message.c_str();
};
