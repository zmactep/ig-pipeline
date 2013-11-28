from unittest import TestCase
from parsers import kabat as kabat


class TestKabat(TestCase):
    def test_empty_sequence(self):
        self.assertEqual(kabat.parse([]), {})

    def test_correct_sequence(self):
        self.assertEqual(kabat.parse(["IGHV4-31*08_IGHD2-15*01{0}_IGHJ4*02{2}	1	90",
                                      "\tIGHV4-34*01_IGHD1-20*01{0}_IGHJ6*02{1}	1	90\t"]),
                         {"IGHV4-31*08_IGHD2-15*01{0}_IGHJ4*02{2}": [(0, 89)],
                          "IGHV4-34*01_IGHD1-20*01{0}_IGHJ6*02{1}": [(0, 89)]})

    def test_malformed_sequence(self):
        self.assertRaises(kabat.ParseException, lambda: kabat.parse(["IGHV4-31*08_IGHD2-15*01{0}_IGHJ4*02{2}	1.0	90.0"]))
        self.assertRaises(kabat.ParseException, lambda: kabat.parse(["IGHV4-31*08_IGHD2-15*01{0}_IGHJ4*02{2}	1 NaN"]))
        self.assertRaises(kabat.ParseException, lambda: kabat.parse(["IGHV4-31*08_IGHD2-15*01{0}_IGHJ4*02{2}	1   90  180"]))
        self.assertRaises(kabat.ParseException, lambda: kabat.parse(["1	90"]))
