#pragma once

#include <algorithm>
#include <limits>
#include <vector>
#include <set>

typedef unsigned char byte;

template<class T>
class kstatistics {
public:
    typedef T link_type;

    kstatistics(std::vector<byte> alphabet, int _k)
        : m_alphabet(alphabet), m_k(_k)
    {
        size_t power = 1;
        for (size_t i = 0; i != m_k; ++i)
        {
            m_pow_cache.push_back(power);
            power *= alphabet.size();
        }
        power *= alphabet.size();
        m_data = std::vector<std::set<link_type>>(power, std::set<link_type>());
    }

    template <class Iterator>
    bool add(Iterator begin, Iterator end, link_type link)
    {
        if (begin > end - m_k)
        {
            return false;
        }
        size_t hash_val = hash(begin, begin + m_k);
        if (hash_val == std::numeric_limits<size_t>::max())
        {
            return false;
        }
        m_data[hash_val].insert(link);
        return true;
    }

    template <class Iterator>
    const std::set<link_type>* get(Iterator begin, Iterator end) const
    {
        if (begin > end - m_k)
        {
            return nullptr;
        }
        size_t hash_val = hash(begin, begin + m_k);
        if (hash_val == std::numeric_limits<size_t>::max())
        {
            return nullptr;
        }
        return &m_data[hash_val];
    }

    template <class Iterator>
    size_t hash(Iterator begin, Iterator end) const
    {
        size_t hash_val = 0;
        if (begin > end - m_k)
        {
            return std::numeric_limits<size_t>::max();
        }
        Iterator iter = begin;
        size_t power = 0;
        while (iter != end) {
            std::vector<byte> ::const_iterator fiter =
                    std::find(m_alphabet.begin(), m_alphabet.end(), (byte)*iter);
            if (fiter == m_alphabet.end())
            {
                return std::numeric_limits<size_t>::max();
            }
            hash_val += m_pow_cache[power] * (fiter - m_alphabet.begin());
            power++;
            ++iter;
        }
        return hash_val;
    }

    bool inAlphabet(byte symbol) const
    {
        if (std::find(m_alphabet.begin(), m_alphabet.end(), symbol) != m_alphabet.end())
        {
            return true;
        }
        return false;
    }

private:
    std::vector<byte>                m_alphabet;
    std::vector<std::set<link_type>> m_data;
    std::vector<size_t>              m_pow_cache;
    int                              m_k;
};
