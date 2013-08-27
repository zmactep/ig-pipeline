#include <ctime>
#include <iostream>
#include <algorithm>
#include "unittest.h"

using namespace std;

typedef unsigned char byte;

void test()
{
    time_t s = clock();

    test_trie();
    test_anno();
    test_kstat();
    test_contig();
    test_contig2();
    test_search();
    test_fasta();
    test_fasta_push();
    test_alicont();
    test_contig_alicont();
    test_contig_alicont2();
    test_contig_alicont3();
    s = clock() - s;
    std::cout << "Test time: " << (double) s / CLOCKS_PER_SEC << std::endl;
}

int main()
{
    test();

    return 0;
}
