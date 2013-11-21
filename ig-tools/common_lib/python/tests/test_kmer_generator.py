from unittest import TestCase

from utils import kmer_generator

class TestKmerGenerator(TestCase):
    def test_get_string_kmers(self):
        self.assertRaises(ValueError, lambda: list(kmer_generator.get_string_kmers("", 0)))
        self.assertRaises(ValueError, lambda: list(kmer_generator.get_string_kmers("abcd", 5)))
        self.assertEqual(list(kmer_generator.get_string_kmers("abcd", 4)), ["abcd"])
        self.assertEqual(list(kmer_generator.get_string_kmers("abcde", 3)), ["abc", "bcd", "cde"])

    def test_get_sequence_kmers(self):
        # this method is also indirectly tested in test_get_string_kmers
        self.assertRaises(ValueError, lambda: kmer_generator.get_sequence_kmers([], 0))
        self.assertEqual(list(kmer_generator.get_sequence_kmers(['a', 'b', 'c'], 2),), [('a', 'b'), ('b', 'c')])