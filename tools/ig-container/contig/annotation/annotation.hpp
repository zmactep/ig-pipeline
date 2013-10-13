#pragma once

#include <string>
#include <vector>
#include <map>

#include "annotation_node.hpp"
#include "annotation_record.hpp"

template <class T, template <class> class Property, class LabelType = std::string>
class annotation
{
public:
    typedef annotation_record<T, Property>    record_type;
    typedef std::pair<LabelType, record_type> value_type;
    typedef typename record_type::data_type   data_type;

    struct index
    {
        size_t record;
        size_t letter;

        index(size_t r = 0, size_t l = 0) : record(r), letter(l)
        {
        }
    };

    annotation(size_t length = 0)
    {
        m_data.reserve(length);
    }

    record_type & add(LabelType const & label, size_t length = 0)
    {
        m_data.push_back(std::make_pair(label, record_type(length)));
        return (*this)[size() - 1];
    }

    record_type & operator[](size_t i)
    {
        return m_data[i].second;
    }

    record_type & operator[](size_t i) const
    {
        return m_data[i].second;
    }

    data_type & operator[](index i)
    {
        return m_data[i.record].second[i.letter];
    }

    const data_type & operator[](index i) const
    {
        return m_data[i.record].second[i.letter];
    }

    LabelType labelOf(size_t i) const
    {
        return m_data[i].first;
    }

    LabelType labelOf(index i) const
    {
        return m_data[i.record].first;
    }

    index last()
    {
        return index(size() - 1, (*this)[size() - 1].size() - 1);
    }

    size_t size() const
    {
        return m_data.size();
    }

    record_type & getRecordByLabel(LabelType const & label)
    {
        for (typename std::vector<value_type>::iterator i = m_data.begin(); i != m_data.end(); ++i)
        {
            if (i->first == label)
            {
                return i->second;
            }
        }
        throw std::range_error("No such name!");
    }

private:
    std::vector<value_type> m_data;
};
