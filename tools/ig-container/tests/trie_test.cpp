#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>

#include "../contig/trie/trie.hpp"

#include "../unittest.h"

template <class T>
void add_test1( trie<T> * t )
{
    std::string s1 = "ACCGT";
    std::string s2 = "ACTT";
    std::string s3 = "CAATGC";
    std::vector<std::string> v;
    v.push_back(s1);
    v.push_back(s2);
    v.push_back(s3);

    for (auto s : v)
    {
        auto iter = t->begin();
        for (auto i : s)
        {
            iter = t->insert(iter, (byte)i);
        }
    }
}

template <class T>
void add_test2( trie<T> * t )
{
    std::string s1 = "ACCGT";
    std::string s2 = "ACTT";
    std::string s3 = "CAATGC";
    std::vector<std::string> v;
    v.push_back(s1);
    v.push_back(s3);
    v.push_back(s2);

    for (auto s : v)
    {
        auto iter = t->begin();
        for (auto i : s)
        {
            iter = t->insert(iter, (byte)i);
        }
    }
}

void iterators_tests()
{
    trie<bool> t;
    trie<bool> s;

    add_test1(&t);
    for (trie<bool>::const_iterator i = t.begin(); i != t.end(); ++i)
    {
        std::cout << i.symbol();
    }
    std::cout << std::endl;

    add_test2(&s);
    for (auto i = s.begin(); i != s.end(); ++i)
    {
        std::cout << i.symbol();
    }
    std::cout << std::endl;

    for (size_t k = 0; k < t.size(); ++k)
    {
        std::cout << (t.begin() + k).symbol();
    }
    std::cout << std::endl;
}

void test_trie()
{
    iterators_tests();
}
