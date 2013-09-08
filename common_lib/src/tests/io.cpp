#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <gtest/gtest.h>
#include "nomenclature.h"

TEST(Model, NomenclatureTest) {
	const std::string line = "A1	1	69	70	117	118	162	163	183	184	279	VK	0";
	Nomenclature a(line);
	//TODO //FIXME
	EXPECT_TRUE(a.getRefName() == "A1");
	EXPECT_TRUE(a.getFR1begin() == 1);
	EXPECT_TRUE(a.getCDR1begin() == 70);
	EXPECT_TRUE(a.getFR2begin() == 118);
	EXPECT_TRUE(a.getCDR2begin() == 163);
	EXPECT_TRUE(a.getFR3begin() == 184);
	EXPECT_TRUE(a.getFR3end() == 279);
}

TEST(IO, SimpleTest) {
	std::ifstream input("./bin/data/sample.txt");
	ASSERT_TRUE(input.is_open() == true);
	std::string line;
	std::getline(input, line);
	ASSERT_TRUE(line == "test");
	input.close();
}

int main(int argc, char **argv) {

	::testing::InitGoogleTest( &argc, argv );
	return RUN_ALL_TESTS();
}
