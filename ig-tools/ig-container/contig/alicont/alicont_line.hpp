#pragma once

#include <limits>
#include <vector>
#include <algorithm>

enum struct path_value
{
    up, left, diag
};

class alicont_line
{
public:
    typedef std::vector<path_value> path_type;

    alicont_line( size_t length, size_t depth )
        : m_score_line(length, path_type(depth, path_value::left)),
          m_score(length, std::numeric_limits<int>::min())
    {
    }

    int operator[](size_t i) const
    {
        return m_score[i];
    }

    const path_type & path(size_t i) const
    {
        return m_score_line[i];
    }

    void set(size_t i, path_type const & prev_path, path_value const & val,
             int score_m, int score_x, int score_y)
    {
        path_type tmp = prev_path;
        tmp.push_back(val);

        m_score_line[i] = tmp;
        m_score[i] = score;
    }

private:
    std::vector<path_type> m_score_line;
    std::vector<int>       m_score;
};
