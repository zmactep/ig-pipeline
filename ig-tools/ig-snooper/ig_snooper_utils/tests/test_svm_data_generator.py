from unittest import TestCase
from itertools import islice
from ig_snooper_utils.svm_data_generator import *


class TestSvmDataGenerator(TestCase):
    def test_get_labels(self):
        self.assertEqual(list(islice(get_labels([]), 3)), [0, 0, 0])
        self.assertEqual(list(get_labels([(0, 1), (2, 5)])), [0, 0, 1, 1, 1, 1])
        self.assertRaises(ValueError, lambda: get_labels([(1, 1)]))
        self.assertRaises(ValueError, lambda: get_labels([(0, 1), (3, 4)]))
        self.assertRaises(ValueError, lambda: get_labels([(0, 2), (1, 3)]))

    def test_get_kmers(self):
        self.assertEqual(list(get_kmers("a", 1, False)), ["a"])
        self.assertEqual(list(get_kmers("a", 1, True)), ["a"])
        self.assertEqual(list(get_kmers("a", 3, True)), [CAP_CHAR + 'a' + CAP_CHAR])

    def test_get_dataset(self):
        self.assertEqual(list(get_dataset('a', [], 3, True)), [(0, CAP_CHAR + 'a' + CAP_CHAR)])
        self.assertEqual(list(get_dataset('ab', [(0, 0), (1, 1)], 3, True)), [(0, CAP_CHAR + 'ab'), (1, 'ab' + CAP_CHAR)])
