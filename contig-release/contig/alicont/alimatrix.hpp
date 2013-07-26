#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <exception>

#include "alialgo.hpp"

// TODO: implement iterators for this container to
//       enable begin-end and c++11 syntax
class alimatrix
{
public:
    alimatrix(const std::vector<simple_matrix2i*>* data)
        : m_data(data)
    {
    }

    std::vector<int> & operator[](size_t index)
    {
        for (std::vector<simple_matrix2i*>::const_iterator i = m_data->begin();
                                                           i != m_data->end();
                                                           ++i)
        {
            if ((*i)->size() < index)
            {
                index -= (*i)->size();
            }
            else
            {
                return (*i)->at(index);
            }
        }

        throw std::out_of_range("Shit happend");
    }

    const std::vector<int> & operator[](size_t index) const
    {
        for (std::vector<simple_matrix2i*>::const_iterator i = m_data->begin();
                                                           i != m_data->end();
                                                           ++i)
        {
            if ((*i)->size() <= index)
            {
                index -= (*i)->size();
            }
            else
            {
                return (*i)->at(index);
            }
        }

        throw std::out_of_range("alimatrix::operator[]");
    }

    size_t size() const
    {
        size_t sum = 0;
        for (std::vector<simple_matrix2i*>::const_iterator i = m_data->begin();
                                                           i != m_data->end();
                                                           ++i)
        {
            sum += (*i)->size();
        }
        return sum;
    }

    void print() const
    {
        for (size_t i = 0; i < size(); ++i)
        {
            const std::vector<int> & v = (*this)[i];
            for (size_t j = 0; j < v.size(); ++j)
            {
                std::cout << std::setw(3) << v[j] << " ";
            }
            std::cout << std::endl;
        }
    }

private:
    const std::vector<simple_matrix2i*>* m_data;
};
