#pragma once

#include <memory>
#include <vector>
#include <map>

#include "trie_node.hpp"
#include "trie_impl.hpp"
#include "trie_cont.hpp"

#include "trie_const_iterator.hpp"
#include "trie_iterator.hpp"

class trie_impl;

template <class T>
class trie
{
public:
    template <class> friend class trie;
    friend class trie_iterator<T>;
    friend class trie_const_iterator<T>;
    typedef trie_node::index_type        index_type;
    typedef trie_const_iterator<T>		 const_iterator;
    typedef trie_iterator<T>             iterator;

    trie() : m_trie(new trie_impl()), m_pure(true)
    {
        m_cont.push(m_trie->root());
    }

    template <class S>
    trie(trie<S> const & parent_trie) : m_trie(parent_trie.m_trie),
                                        m_cont(parent_trie.m_cont),
                                        m_pure(false)
    {
    }

    ~trie()
    {
        if (m_pure)
        {
            delete m_trie;
        }
    }

    iterator insert(iterator i, byte symbol)
    {
        trie_node* new_node = m_trie->insert(m_cont.nodeOf(i.index()), symbol);
        if (new_node->id() == m_cont.size())
        {
            return iterator(this, m_cont.push(new_node));
        }
        return iterator(this, new_node->id());
    }

    iterator begin()
    {
        return iterator(this, 0);
    }

    iterator end()
    {
        return iterator(this, -1);
    }

    const_iterator begin() const
    {
        return const_iterator(this, 0);
    }

    const_iterator end() const
    {
        return const_iterator(this, -1);
    }

    size_t size() const
    {
        return m_cont.size();
    }

    iterator iter(size_t i)
    {
        return iterator(this, i);
    }

    const_iterator iter(size_t i) const
    {
        return const_iterator(this, i);
    }

private:
    T & operator[](size_t i)
    {
        return m_cont[i];
    }

    T & operator[](size_t i) const
    {
        return m_cont[i];
    }

    trie_node* get(size_t i) const
    {
        return m_cont.nodeOf(i);
    }

    size_t next(size_t i, byte symbol) const
    {
        trie_node* node = m_cont.nodeOf(i)->next(symbol);
        if (node != nullptr)
        {
            return node->id();
        }

        return -1;
    }

    size_t prev(size_t i) const
    {
        trie_node* node = m_cont.nodeOf(i)->prev();
        if (node != nullptr)
        {
            return node->id();
        }

        return -1;
    }

    size_t dfsNext(size_t i) const
    {
        return m_trie->dfs_cache[i];
    }

    size_t dfs(size_t id, size_t last)
    {
        m_trie->dfs_cache[last] = id;
        last = id;
        trie_node::children_t* children = m_cont.nodeOf(id)->getChildren();
        for (trie_node::children_t::iterator iter = children->begin(); iter != children->end(); ++iter)
        {
            last = dfs(iter->second->id(), last);
        }

        return last;
    }

    void createDfsCache()
    {
        if (m_trie->dfs_cache.size() == m_cont.size())
            return;

        m_trie->dfs_cache = std::vector<size_t>(m_cont.size(), -1);
        dfs(0, 0);
    }

private:
    trie_impl*	         m_trie;
    trie_cont<T>         m_cont;
    //std::vector<size_t>  m_dfs_cache;
    bool                 m_pure;
};
