#pragma once

#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "contig/annotation/annotation.hpp"
#include "contig/trie/trie.hpp"
#include "contig/kstat/kstat.hpp"
#include "contig/contig.hpp"

#include "utils/fasta_reader.h"

typedef unsigned char byte;

struct Alphabet
{
    static std::vector<byte> getAlphabet()
    {
        std::vector<byte> v;
        v.push_back('A');
        v.push_back('C');
        v.push_back('G');
        v.push_back('T');
        v.push_back('N');
        return v;
        //return std::vector<byte> {'A', 'C', 'G', 'T', 'N'};
    }
};

template <class T>
struct RegionProp
{
    size_t region_id;
    std::string name;

    RegionProp(size_t id = 0) : region_id(id)
    {
    }
};

std::pair<Read, size_t> import_data(std::string const & path,
                                    contig<Alphabet, RegionProp> & my_contig);

void test_trie();
void test_kstat();
void test_anno();
void test_contig();
void test_contig2();
void test_search();
void test_fasta();
void test_fasta_push();
void test_alicont();
void test_contig_alicont();
void test_contig_alicont2();
