#pragma once

#include <vector>
#include <map>
#include <string>
#include <limits>

#include "score_matrix.hpp"

#define max3(a,b,c) std::max(a, std::max(b, c))

typedef unsigned char byte;

typedef std::vector<std::vector<int>> simple_matrix2i;

struct alignment_result
{
    std::string first;
    std::string second;
    int         score;

    size_t      target_id;

    alignment_result() : score(0)
    {
    }

    alignment_result(std::string const & f,
                     std::string const & s,
                     int sc)
        : first(f), second(s), score(sc), target_id(std::numeric_limits<size_t>::max())
    {
    }

    bool operator<(alignment_result const & ar) const
    {
        return score < ar.score;
    }
};

template <class StringOrRope>
simple_matrix2i needleman_wunsch(StringOrRope const & seq1, StringOrRope const & seq2,
                                 int gap, score_matrix const & matrix,
                                 simple_matrix2i* prev_matrix = nullptr)
{
    std::vector<int>* prev;
    simple_matrix2i   d;
    if (prev_matrix == nullptr)
    {
        d.push_back(std::vector<int>(seq1.size() + 1, 0));
        for (size_t i = 0; i < d[0].size(); ++i)
        {
            d[0][i] = i * gap;
        }
        prev = &d[0];
    }
    else
    {
        prev = &(prev_matrix->back());
    }

    for (size_t i = 1; i < seq2.size() + 1; ++i)
    {
        byte t = seq2[i - 1];
        std::vector<int> current(seq1.size() + 1, 0);
        current[0] = prev->at(0) + gap;
        for (size_t j = 0; j < seq1.size(); ++j)
        {
            byte s = seq1[j];
            current[j + 1] = max3(prev->at(j) + matrix(s, t),
                                  prev->at(j + 1) + gap,
                                  current[j] + gap);
        }
        d.push_back(current);
        prev = &(d.back());
    }

    return d;
}

template <class StringOrRope, class Container2d>
alignment_result traceback(StringOrRope const & seq1, StringOrRope const & seq2,
                           Container2d const & score, int gap,
                           score_matrix const & matrix)
{
    size_t i = seq2.size();
    size_t j = seq1.size();

    std::string res1;
    std::string res2;

    while (i != 0 || j != 0)
    {
        byte s = seq2[i - 1];
        byte t = seq1[j - 1];
        int sc = score[i][j];

        if (i != 0 && sc == score[i - 1][j] + gap)
        {
            i -= 1;
            res1.push_back('-');
            res2.push_back(s);
        }
        else if (j != 0 && sc == score[i][j - 1] + gap)
        {
            j -= 1;
            res1.push_back(t);
            res2.push_back('-');
        }
        else if (sc == score[i - 1][j - 1] + matrix(s, t))
        {
            i -= 1;
            j -= 1;
            res1.push_back(t);
            res2.push_back(s);
        }
    }

    std::reverse(res1.begin(), res1.end());
    std::reverse(res2.begin(), res2.end());
    return alignment_result(res1, res2, score[seq2.size()][seq1.size()]);
}
