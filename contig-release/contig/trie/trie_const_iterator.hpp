#pragma once

#include <iterator>

#include "trie.hpp"

template <class> class trie;

template <class T>
class trie_const_iterator : public std::iterator<std::bidirectional_iterator_tag,
                                                 T, std::ptrdiff_t,
                                                 const T*, const T&>
{
public:
    friend class trie<T>;

    trie_const_iterator() : m_trie(nullptr), m_current(-1)
    {
    }

    trie_const_iterator(trie<T>* t, size_t current = 0) : m_trie(t), m_current(current)
    {
    }

    trie_const_iterator next(byte symbol) const
    {
        return trie_const_iterator(m_trie, m_trie->next(m_current, symbol));
    }

    trie_const_iterator prev() const
    {
        return trie_const_iterator(this, m_trie->prev(m_current));
    }

    trie_const_iterator & operator++()
    {
        m_trie->createDfsCache();
        m_current = m_trie->dfsNext(m_current);
        return *this;
    }

    trie_const_iterator operator++(int)
    {
        m_trie->createDfsCache();
        trie_const_iterator tmp = *this;
        m_current = m_trie->dfsNext(m_current);
        return tmp;
    }

    trie_const_iterator & operator--()
    {
        m_current = m_trie->prev(m_current);
        return *this;
    }

    trie_const_iterator operator--(int)
    {
        trie_const_iterator tmp = *this;
        m_current = m_trie->prev(m_current);
        return tmp;
    }

    trie_const_iterator operator+(size_t i) const
    {
        m_trie->createDfsCache();
        size_t new_current = m_current;
        for (size_t k = 0; k < i; ++k)
        {
            new_current = m_trie->dfsNext(new_current);
        }
        return trie_const_iterator(m_trie, new_current);
    }

    trie_const_iterator operator-(size_t i) const
    {
        size_t new_current = m_current;
        for (size_t k = 0; k < i; ++k)
        {
            new_current = m_trie->prev(new_current);
        }
        return trie_const_iterator(m_trie, new_current);
    }

    T & operator*() const
    {
        return (*m_trie)[m_current];
    }

    T* operator->() const
    {
        return &(**this);
    }

    bool operator==(const trie_const_iterator<T> iter) const
    {
        if (m_current < m_trie->size())
            return m_current == iter.m_current;
        return iter.m_current > m_trie->size();
    }

    bool operator!=(const trie_const_iterator<T> iter) const
    {
        return !(*this == iter);
    }

    byte symbol() const
    {
        return m_trie->get(m_current)->symbol();
    }

    size_t index() const
    {
        return m_current;
    }

    bool leaf() const
    {
        return m_trie->get(m_current)->getChildren()->size() == 0;
    }

    bool fork() const
    {
        return m_trie->get(m_current)->getChildren()->size() > 1;
    }

protected:
    trie<T>*  m_trie;
    size_t    m_current;
};
