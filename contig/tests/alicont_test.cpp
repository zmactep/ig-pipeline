#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <ctime>

#include "../contig/contig.hpp"
#include "../contig/alicont/alicont.hpp"

#include "../unittest.h"

std::pair<Read, size_t> import_data_a(std::string const & path,
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
        std::vector<RegionProp<Alphabet>> v(tmp_read.seq.size());
        for (size_t i = 0; i < v.size(); ++i)
        {
            v[i].region_id = rand() % 7;
        }
        std::cout << std::endl;
        my_contig.push(tmp_read.seq.begin(), tmp_read.seq.end(), v.begin(), tmp_read.name);
        real_length += tmp_read.seq.size();
        count++;
    }
    std::cout << "Import ended (" << count << ")" << std::endl;

    return std::make_pair(tmp_read, real_length);
}

void test_alicont()
{
    score_matrix m("../data/BLOSUM62.txt");
    alicont a("MEANLY", -5, std::move(m));
    simple_matrix2i m1 = std::move(a.score("PLEA"));
    a.push("PLEA", &m1);
    simple_matrix2i m2 = std::move(a.score("SANT"));
    a.push("SANT", &m2);
    simple_matrix2i m3 = std::move(a.score("LY"));
    a.push("LY", &m3);

    std::cout << a.target() << std::endl;
    a.print();
    alignment_result res = a.alignment();
    std::cout << res.first << std::endl << res.second << std::endl;

    a.pop();
    m3 = std::move(a.score("PYR"));
    a.push("PYR", &m3);

    std::cout << a.target() << std::endl;
    a.print();
    res = a.alignment();
    std::cout << res.first << std::endl << res.second << std::endl;
}

void test_contig_alicont()
{
    contig<Alphabet, RegionProp> my_contig("CONTIG-TEST", Alphabet::getAlphabet());
    std::string s1 = "ACGTAG";
    std::string s2 = "ACGTCGT";
    std::string s3 = "ACAAT";
    std::string s4 = "ACAAC";
    std::string s5 = "CCGTTG";

    my_contig.push(s1.begin(), s1.end(), "s1");
    my_contig.push(s2.begin(), s2.end(), "s2");
    my_contig.push(s3.begin(), s3.end(), "s3");
    my_contig.push(s4.begin(), s4.end(), "s4");
    my_contig.push(s5.begin(), s5.end(), "s5");

    std::string query = "GCGTTG";
    score_matrix m("../data/BLOSUM62.txt");
    auto r = my_contig.align(query.begin(), query.end(), -5, m, 3);
    for (auto i : r)
    {
        std::cout << i.score << std::endl << i.first
                             << std::endl << i.second
                             << std::endl << std::endl;
    }
    std::cout << std::endl;
}

void test_contig_alicont2()
{
    contig<Alphabet, RegionProp> my_contig("CONTIG-TEST", Alphabet::getAlphabet());
    //import_data("../data/germline/human/VH.fasta", my_contig);
    import_data("../data/other/VH_corrected.fasta", my_contig);
    std::string query = "CAGGTTCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTTTCCTGCAAGGCTTCTGGATACACCTTCACTAGCTATGCTATGCATTGGGTGCGCCAGGCCCCCGGACAAAGGCTTGAGTGGATGGGATGGAGCAACGCTGGCAATGGTAACACAAAATATTCACAGGAGTTCCAGGGCAGAGTCACCATTACCAGGGACACATCCGCGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACATGGCTGTGTATTACTGTGCGAGAGA";

    score_matrix m("../data/NUC4.4.txt");
    time_t t = clock();
    auto r = my_contig.align_sc(query.begin(), query.end(), -5, m, 0.8);
    std::cout << "Alignment time: " << (double)(clock() - t) / CLOCKS_PER_SEC << std::endl;
    for (auto i : r)
    {
        std::cout << i.score  << " (" << i.target_id << ")" << std::endl
                  << i.first  << std::endl
                  << i.second << std::endl << std::endl;
    }
    std::cout << std::endl;
}

void test_contig_alicont3()
{
    contig<Alphabet, RegionProp> my_contig("CONTIG-TEST", Alphabet::getAlphabet());
    import_data_a("../data/germline/human/VH.fasta", my_contig);
    std::string query = "CAGGTTTCTGGGGCTGAGGTGAAGAAGCCTGGGGGCAAGGCTTCTGGATACACCTTCACTAGCTATGCTATGCATTGGGTGCGCCAGGCCCCCGGACAAAGGCTTGAGTGGATGGGATGGAGCAACGCTGGCAATGGTAACACAAAATATTCACAGGAGTTCCAGGGCAGAGTCACCATTACCAGGGACACATCCGCGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACATGGCTGTGTATTACTGTGCGAGAGA";

    score_matrix m("../data/NUC4.4.txt");
    time_t t = clock();
    auto r = my_contig.annotate(query.begin(), query.end(), -5, m);
    std::cout << "Alignment time: " << (double)(clock() - t) / CLOCKS_PER_SEC << std::endl;
    for (size_t i = 0; i < r.size(); ++i)
    {
        std::cout << r[i].region_id;
    }
    std::cout << std::endl;
}
