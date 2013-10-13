#include <iostream>
#include <vector>
#include <string>
#include <gtest/gtest.h>
#include "tokenizer.h"

TEST(Utils, TokenizerTest) {
	const std::string line = "Hello World";
	std::vector<std::string> tokens = Tokenizer::tokenize(line);
	ASSERT_TRUE(tokens.size() == 2);
	ASSERT_TRUE(tokens[0] == "Hello");
	ASSERT_TRUE(tokens[1] == "World");
}

int main(int argc, char **argv) {

	::testing::InitGoogleTest( &argc, argv );
	return RUN_ALL_TESTS();
}
