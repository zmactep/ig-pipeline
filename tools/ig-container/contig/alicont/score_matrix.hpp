#pragma once

#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>

class score_matrix
{
public:
    score_matrix()
    {
    }

    score_matrix(score_matrix const & sm)
        : m_map(sm.m_map)
    {
    }

    score_matrix(score_matrix && sm)
        : m_map(std::move(sm.m_map))
    {
    }

    score_matrix(std::string const & filename)
    {
        load_matrix(filename);
    }

    ~score_matrix()
    {
    }

    const int & operator()(char _x, char _y) const
    {
        return m_map.at(std::make_pair(_x, _y));
    }

    void load_matrix(std::string const & filename)
    {
        m_map.clear();

        std::ifstream infile(filename);
        std::string line;
        std::string amino;
        char current_letter_first;
        char current_letter_second;
        int value;
        int x=0;
        int y=0;

        while (std::getline(infile, line))
        {
            if (line[0] == '#')
            {
                continue;
            }
            if(line[0] == ' ')
            {
                amino = line;
                amino.erase(std::remove_if(amino.begin(), amino.end(), isspace), amino.end());
                continue;
            }

            x = 0;
            std::stringstream converter(&line[0]); // &line[1] to skip leading letter

            converter >> current_letter_first;

            while (converter >> value) {
                current_letter_second = amino[x];
                std::pair<char, char> key = std::make_pair(current_letter_first, current_letter_second);
                m_map.insert(std::make_pair(key, value));
                x++;
            }
            ++y;
        }
    }

private:
    std::map<std::pair<char, char>, int> m_map;
};
