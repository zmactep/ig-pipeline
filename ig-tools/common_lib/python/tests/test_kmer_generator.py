from unittest import TestCase

from utils import kmer_generator

class TestKmerGenerator(TestCase):
    def test_empty(self):
        self.assertRaises(ValueError, lambda: list(kmer_generator.get_string_kmers("", 0)))

    def test_shorter_than_window(self):
        self.assertRaises(ValueError, lambda: list(kmer_generator.get_string_kmers("abcd", 5)))

    def test_one(self):
        self.assertEqual(list(kmer_generator.get_string_kmers("abcd", 4)), ["abcd"])

    def test_several(self):
        self.assertEqual(list(kmer_generator.get_string_kmers("abcde", 3)), ["abc", "bcd", "cde"])