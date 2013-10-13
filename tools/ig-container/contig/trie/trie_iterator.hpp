#pragma once

#include <iterator>

#include "trie.hpp"
#include "trie_const_iterator.hpp"

template <class> class trie;

template <class T>
class trie_iterator : public trie_const_iterator<T>
{
public:
    friend class trie<T>;

    trie_iterator() : trie_const_iterator<T>()
    {
    }

    trie_iterator(trie<T>* t, size_t current = 0)
        : trie_const_iterator<T>(t, current)
    {
    }

    trie_iterator(trie_const_iterator<T> t)
    {
        this->m_trie = t.m_trie;
        this->m_current = t.m_current;
    }

    trie_const_iterator<T> const_iter() const
    {
        return trie_const_iterator<T>(this->m_trie, this->m_current);
    }

    trie_iterator next(byte symbol) const
    {
        return trie_iterator(this->m_trie,
                             this->m_trie->next(this->m_current, symbol));
    }

    trie_iterator prev() const
    {
        return trie_iterator(this->m_trie, this->m_trie->prev(this->m_current));
    }

    trie_iterator & operator++()
    {
        this->m_trie->createDfsCache();
        this->m_current = this->m_trie->dfsNext(this->m_current);
        return *this;
    }

    trie_iterator operator++(int)
    {
        this->m_trie->createDfsCache();
        trie_iterator tmp = *this;
        this->m_current = this->m_trie->dfsNext(this->m_current);
        return tmp;
    }

    trie_iterator & operator--()
    {
        this->m_current = this->m_trie->prev(this->m_current);
        return *this;
    }

    trie_iterator operator--( int)
    {
        trie_iterator tmp = *this;
        this->m_current = this->m_trie->prev(this->m_current);
        return tmp;
    }

    trie_iterator operator+(size_t i) const
    {
        this->m_trie->createDfsCache();
        size_t new_current = this->m_current;
        for (size_t k = 0; k < i; ++k)
        {
            new_current = this->m_trie->dfsNext(new_current);
        }
        return trie_iterator(this->m_trie, new_current);
    }

    trie_iterator operator-(size_t i) const
    {
        size_t new_current = this->m_current;
        for (size_t k = 0; k < i; ++k)
        {
            new_current = this->m_trie->prev(new_current);
        }
        return trie_iterator(this->m_trie, new_current);
    }

    T & operator*() const
    {
        return (*(this->m_trie))[this->m_current];
    }

    T* operator->() const
    {
        return &(**this);
    }

    bool operator==(trie_iterator<T> iter) const
    {
        if (this->m_current < this->m_trie->size())
            return this->m_current == iter.m_current;
        return iter.m_current > this->m_trie->size();
    }

    bool operator!=(trie_iterator<T> iter) const
    {
        return !(*this == iter);
    }
};
