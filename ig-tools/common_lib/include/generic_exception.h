#ifndef GENERICEXCEPTION_H_
#define GENERICEXCEPTION_H_

#include <exception>
#include <string>

class GenericException : public std::exception {
public:
	GenericException(const std::string & m) {
		this->message = m;
	}

	~GenericException() throw (){};

	const char* what() const throw (){
		return message.c_str();
	};

private:
	std::string message;
};

#endif /* GENERICEXCEPTION_H_ */
