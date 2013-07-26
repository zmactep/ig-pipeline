#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "../contig/annotation/annotation.hpp"
#include "../contig/trie/trie.hpp"
#include "../contig/kstat/kstat.hpp"
#include "../contig/contig.hpp"

#include "../utils/fasta_reader.h"

#include "../unittest.h"

std::pair<Read, size_t> import_data(std::string const & path,
                                    contig<Alphabet, RegionProp> & my_contig)
{
    FastaReader fr(path);
    Read tmp_read;
    std::cout << "Import started" << std::endl;
    size_t count = 0;
    size_t real_length = 0;
    while (!fr.is_eof())
    {
        fr >> tmp_read;
        my_contig.push(tmp_read.seq.begin(), tmp_read.seq.end(), tmp_read.name);
        real_length += tmp_read.seq.size();
        count++;
    }
    std::cout << "Import ended (" << count << ")" << std::endl;

    return std::make_pair(tmp_read, real_length);
}

void test_fasta()
{
    time_t c = clock();
    std::pair<Read, size_t> p;
    size_t real_length = 0;
    contig<Alphabet, RegionProp> my_contig("CONTIG-TEST", Alphabet::getAlphabet());
//    p = import_data("../../../../data/germline/human/VH.fasta", my_contig);
//    real_length += p.second;
    p = import_data("/home/mactep/Data/NGS-llama/out2/VL/VL.fasta", my_contig);
    real_length += p.second;
//    p = import_data("/home/mactep/Data/NGS-llama/out2/VK/VK.fasta", my_contig);
//    real_length += p.second;

    std::cout << "FOUND: "
              << my_contig.find(p.first.seq.begin() + 12,
                                p.first.seq.begin() + 22).size()
              << std::endl;
    c = clock() - c;
    std::cout << "Compression: " << (double)my_contig.size() / real_length << std::endl;
    std::cout << "Real time: " << (double)c / CLOCKS_PER_SEC << std::endl;
}

void test_fasta_push()
{
    contig<Alphabet, RegionProp> my_contig("CONTIG-TEST", Alphabet::getAlphabet());
    FastaReader fr("/home/mactep/Data/NGS-llama/out2/VH/VH_corrected.fasta");
    Read tmp_read;
    std::vector<Read> tmp_vec;
    time_t c = clock();
    while (!fr.is_eof())
    {
        fr >> tmp_read;
        tmp_vec.push_back(tmp_read);
    }
    std::cout << "KSeq read: " << (double)(clock() - c) / CLOCKS_PER_SEC << std::endl;

    c = clock();
    for (auto r : tmp_vec)
    {
        my_contig.push(r.seq.begin(), r.seq.end(), r.name);
    }
    std::cout << "Container push: " << (double)(clock() - c) / CLOCKS_PER_SEC << std::endl;
}
