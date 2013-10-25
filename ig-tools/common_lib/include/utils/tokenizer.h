#ifndef TOKENIZER_H_
#define TOKENIZER_H_

#include <vector>
#include <string>
class Tokenizer {
public:
	static std::vector<std::string> tokenize(const std::string& line);
};


#endif /* TOKENIZER_H_ */
