#pragma once

#include <iterator>
#include <vector>
#include "contig.hpp"


template <class T, template <class> class Property, class LabelType>
class contig;

template <class T, template <class> class Property, class LabelType=std::string>
class contig_const_iterator : public std::iterator<std::bidirectional_iterator_tag,
                                                   std::vector<Property<T>>, std::ptrdiff_t,
                                                   const std::vector<Property<T>>*,
                                                   const std::vector<Property<T>>&>
{
public:
    friend class contig<T, Property, LabelType>;
    typedef contig<T, Property, LabelType> contig_t;
    typedef Property<T>                    data_type;

    contig_const_iterator() : m_contig(nullptr), m_trie_iter()
    {
    }

    contig_const_iterator(contig_t* contig, size_t current)
        : m_contig(contig), m_trie_iter(&(m_contig->m_trie), current)
    {
    }

    contig_const_iterator(contig_t* contig, typename contig_t::const_iterator trie_iter)
        : m_contig(contig), m_trie_iter(trie_iter)
    {
    }

    contig_const_iterator(contig_t* contig, typename contig_t::iterator trie_iter)
        : m_contig(contig), m_trie_iter(trie_iter)
    {
    }

    contig_const_iterator next(byte symbol) const
    {
        return contig_const_iterator(m_contig, m_trie_iter.next(symbol));
    }

    contig_const_iterator prev() const
    {
        return contig_const_iterator(m_contig, m_trie_iter.prev());
    }

    contig_const_iterator & operator++()
    {
        ++m_trie_iter;
        return *this;
    }

    contig_const_iterator operator++(int)
    {
        contig_const_iterator tmp = *this;
        ++m_trie_iter;
        return tmp;
    }

    contig_const_iterator & operator--()
    {
        --m_trie_iter;
        return *this;
    }

    contig_const_iterator operator--(int)
    {
        contig_const_iterator tmp = *this;
        --m_trie_iter;
        return tmp;
    }

    contig_const_iterator operator+(size_t i) const
    {
        return contig_const_iterator(m_contig, m_trie_iter + i);
    }

    contig_const_iterator operator-(size_t i) const
    {
        return contig_const_iterator(m_contig, m_trie_iter - i);
    }

    std::vector<data_type> & operator*() const
    {
        return m_contig->getAnnotations(m_trie_iter);
    }

    std::vector<data_type>* operator->() const
    {
        return &(**this);
    }

    bool operator==(const contig_const_iterator<T, Property, LabelType> iter) const
    {
        return m_trie_iter == iter.m_trie_iter;
    }

    bool operator!=(const contig_const_iterator<T, Property, LabelType> iter) const
    {
        return !(*this == iter);
    }

    byte symbol() const
    {
        return m_trie_iter.symbol();
    }

    size_t index() const
    {
        return m_trie_iter.index();
    }

    bool leaf() const
    {
        return m_trie_iter.leaf();
    }

    bool fork() const
    {
        return m_trie_iter.fork();
    }

    std::vector<LabelType> getLabels() const
    {
        return m_contig->getLabels(m_trie_iter);
    }

private:
    contig_t*                   m_contig;
    typename contig_t::iterator m_trie_iter;
};
