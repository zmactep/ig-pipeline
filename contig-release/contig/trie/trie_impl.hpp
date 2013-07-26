#pragma once

#include <vector>

#include "trie_node.hpp"

class trie_impl
{
public:
    trie_impl() : m_counter(0), m_root(new trie_node(nullptr, 0, m_counter++))
    {
    }

    ~trie_impl()
    {
        delete m_root;
    }

    trie_node* insert(trie_node* index, byte value)
    {
        trie_node* new_node = index->push(value, m_counter);
        if (new_node->id() == m_counter)
        {
            m_counter++;
        }

        return new_node;
    }

    size_t size() const
    {
        return m_counter;
    }

    trie_node* root()
    {
        return m_root;
    }

private:
    trie_impl(trie_impl const & t);
    trie_impl & operator=(trie_impl const & t);

private:
    size_t               m_counter;
    trie_node*           m_root;

public:
    std::vector<size_t>  dfs_cache;
};
