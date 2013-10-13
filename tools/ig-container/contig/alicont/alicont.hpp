#pragma once

#include <vector>
#include <string>
#include <map>

#include "score_matrix.hpp"

#include "alialgo.hpp"
#include "alimatrix.hpp"

class alicont
{
public:
    typedef std::pair<std::string, simple_matrix2i*> pair_data;

    alicont(std::string const & line, int gap, score_matrix const & matrix)
        : m_seqs(), m_data(), m_line(line), m_gap(gap), m_matrix(matrix)
    {
    }

    alicont(std::string const & line, int gap, score_matrix && matrix)
        : m_seqs(), m_data(), m_line(line), m_gap(gap), m_matrix(std::move(matrix))
    {
    }

    // To use it with good perfomance use std::move
    simple_matrix2i score(std::string const & seq)
    {
        return needleman_wunsch(m_line, seq, m_gap, m_matrix, top().second);
    }

    // std::move for first argument here would be a good idea too
    void push(std::string const & seq, simple_matrix2i * matrix)
    {
        m_seqs.push_back(seq);
        m_data.push_back(matrix);
    }

    inline pair_data pop()
    {
        if (!m_seqs.size())
        {
            return std::make_pair(std::string(), nullptr);
        }
        pair_data result = std::make_pair(m_seqs.back(), m_data.back());
        m_seqs.pop_back();
        m_data.pop_back();
        return result;
    }

    inline pair_data top()
    {
        if (!m_seqs.size())
        {
            return std::make_pair(std::string(), nullptr);
        }
        pair_data result = std::make_pair(m_seqs.back(), m_data.back());
        return result;
    }

    alignment_result alignment()
    {
        return traceback(m_line, target(), alimatrix(&m_data), m_gap, m_matrix);
    }

    inline int score()
    {
        simple_matrix2i* top = m_data.back();
        size_t i = top->size() - 1;
        size_t j = (*top)[i].size() - 1;

        return (*top)[i][j];
    }

    std::string target()
    {
        std::string s;
        for (std::vector<std::string>::iterator i = m_seqs.begin();
                                                i != m_seqs.end();
                                                ++i)
        {
            s += *i;
        }

        return s;
    }

    void clear()
    {
        m_seqs.clear();
        m_data.clear();
    }

    const std::string & query() const
    {
        return m_line;
    }

    void print() const
    {
        alimatrix(&m_data).print();
    }

private:
    std::vector<std::string>      m_seqs;
    std::vector<simple_matrix2i*> m_data;
    std::string                   m_line;
    int                           m_gap;
    score_matrix                  m_matrix;
};
