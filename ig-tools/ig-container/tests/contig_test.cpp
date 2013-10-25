#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>

#include "../contig/annotation/annotation.hpp"
#include "../contig/trie/trie.hpp"
#include "../contig/kstat/kstat.hpp"
#include "../contig/contig.hpp"

#include "../unittest.h"

void test_search()
{
    contig<Alphabet, RegionProp> my_contig("CONTIG-TEST", Alphabet::getAlphabet());
    std::string s1 = "ACCCGTCGTGCAGCATGCATGCGACTACGCGCA";
    std::string s2 = "ACCCATCGATCTGCGACTACGCGCA";
    std::string s3 = "CGAACGCTCAGCATGCATGCGACTACGCGCA";
    std::string s4 = "ACCCATCGATCTGCGACTACGTTC";
    std::string s5 = "ACTCATCGATCTGCGACTACGAAA";

    my_contig.push(s1.begin(), s1.end(), "s1");
    my_contig.push(s2.begin(), s2.end(), "s2");
    my_contig.push(s3.begin(), s3.end(), "s3");
    my_contig.push(s4.begin(), s4.end(), "s4");
    my_contig.push(s5.begin(), s5.end(), "s5");

    std::string sp = "GATCTGCGACTACG";
    auto result = my_contig.find(sp.begin(), sp.end());
    for(auto i : result)
    {
        contig<Alphabet, RegionProp>::const_aiterator iter = my_contig.aiter(i);
        for (auto c = sp.begin(); c != sp.end() - 1; )
        {
            std::cout << *c << ": ";
            for (auto j : iter.getLabels())
            {
                std::cout << j << " ";
            }
            std::cout << std::endl;
            ++c;
            iter = iter.next(*c);
        }
    }
    std::cout << std::endl;
}

void test_contig2()
{
    contig<Alphabet, RegionProp> my_contig("CONTIG-TEST", Alphabet::getAlphabet());
    std::string s1 = "ACCCGTCGTGCAGCATGCATGCGACTACGCGCA";
    std::string s2 = "ACCCATCGATCTGCGACTACGCGCA";
    std::string s3 = "CGAACGCTCAGCATGCATGCGACTACGCGCA";

    my_contig.push(s1.begin(), s1.end(), "s1");
    my_contig.push(s2.begin(), s2.end(), "s2");
    my_contig.push(s3.begin(), s3.end(), "s3");

    std::cout << "Labels of nodes in DFS:" << std::endl;
    contig<Alphabet, RegionProp>::iterator iter = my_contig.begin();
    while (iter != my_contig.end())
    {
        auto labels = my_contig.getLabels(iter++);
        std::cout << "( ";
        for(auto i : labels)
        {
            std::cout << i << " ";
        }
        std::cout << ")   ";
    }
    std::cout << std::endl;
}

void test_contig()
{
    contig<Alphabet, RegionProp> my_contig("CONTIG-TEST", Alphabet::getAlphabet());
    std::string s1 = "ACCCGTCGTGCAGCATGCATGCGACTACGCGCA";

    contig<Alphabet, RegionProp>::const_iterator iter =
            my_contig.push(s1.begin(), s1.end(), "s1");

    std::cout << "get labels" << std::endl;
    auto labels = my_contig.getLabels(iter);
    std::cout << "( ";
    for(auto i : labels)
    {
        std::cout << i << " ";
    }
    std::cout << ")   " << std::endl;

    while (iter != my_contig.end())
    {
        auto annotations = my_contig.getAnnotations(iter++);
        for (auto i : annotations)
        {
            std::cout << i.region_id << " ";
        }
    }
    std::cout << std::endl;

    std::cout << "with root" << std::endl;
    for (auto i = my_contig.begin(); i != my_contig.end(); ++i)
    {
        auto annotations = my_contig.getAnnotations(i);
        for (auto i : annotations)
        {
            std::cout << i.region_id << " ";
        }
    }
    std::cout << std::endl;

    std::cout << "isomorphic tree creation" << std::endl;
    trie<bool>* my_trie = nullptr;
    my_contig.copyTrie(&my_trie);
    for (auto i : *my_trie)
    {
        std::cout << i << " ";
    }
    delete my_trie;
    std::cout << std::endl;
}
