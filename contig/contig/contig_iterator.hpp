#pragma once

#include <iterator>
#include <string>

#include "contig.hpp"
#include "contig_const_iterator.hpp"

template <class T, template <class> class Property, class LabelType>
class contig;

template <class T, template <class> class Property, class LabelType=std::string>
class contig_iterator : public contig_const_iterator<T, Property, LabelType>
{
public:
    friend class contig<T, Property, LabelType>;
    typedef contig<T, Property, LabelType> contig_t;
    typedef Property<T>                    data_type;

    contig_iterator() : contig_const_iterator<T, Property, LabelType>()
    {
    }

    contig_iterator(contig_t* contig, size_t current = 0)
        : contig_const_iterator<T, Property, LabelType>(contig, current)
    {
    }

    contig_iterator(contig_const_iterator<T, Property, LabelType> ci)
    {
        this->m_contig = ci.m_contig;
        this->m_trie_iter = ci.m_trie_iter;
    }

    contig_const_iterator<T, Property, LabelType> const_iter() const
    {
        return contig_const_iterator<T, Property, LabelType>(this->m_contig,
                                                             this->m_trie_iter);
    }

    contig_iterator next(byte symbol) const
    {
        return contig_iterator(this->m_contig, this->m_trie_iter.next(symbol));
    }

    contig_iterator prev() const
    {
        return contig_iterator(this->m_contig, this->m_trie_iter.prev());
    }

    contig_iterator & operator++()
    {
        ++(this->m_trie_iter);
        return *this;
    }

    contig_iterator operator++(int)
    {
        contig_iterator tmp = *this;
        ++(this->m_trie_iter);
        return tmp;
    }

    contig_iterator & operator--()
    {
        --(this->m_trie_iter);
        return *this;
    }

    contig_iterator operator--(int)
    {
        contig_iterator tmp = *this;
        --(this->m_trie_iter);
        return tmp;
    }

    contig_iterator operator+(size_t i) const
    {
        return contig_iterator(this->m_contig, this->m_trie_iter + i);
    }

    contig_iterator operator-(size_t i) const
    {
        return contig_iterator(this->m_contig, this->m_trie_iter - i);
    }

    std::vector<data_type> & operator*()
    {
        return this->m_contig->getAnnotations(this->m_trie_iter);
    }

    std::vector<data_type>* operator->()
    {
        return &(**this);
    }

    bool operator==(const contig_iterator<T, Property, LabelType> iter) const
    {
        return this->m_trie_iter == iter.m_trie_iter;
    }

    bool operator!=(const contig_iterator<T, Property, LabelType> iter) const
    {
        return !(*this == iter);
    }
};
