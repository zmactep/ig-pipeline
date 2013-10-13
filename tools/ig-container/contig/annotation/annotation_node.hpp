#pragma once

typedef unsigned char byte;

template <class T, template <class> class Property>
class annotation_node
{
public:
    typedef Property<T> data_type;

    annotation_node(size_t index, byte symbol)
        : m_data(Property<T>()), m_trie_index(index), m_symbol(symbol)
    {
    }

    data_type & operator*()
    {
        return m_data;
    }

    const data_type & operator*() const
    {
        return m_data;
    }

    size_t index() const
    {
        return m_trie_index;
    }

    byte symbol() const
    {
        return m_symbol;
    }

private:
    Property<T> m_data;
    size_t      m_trie_index;
    byte        m_symbol;
};
