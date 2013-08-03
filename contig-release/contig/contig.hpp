#pragma once

#include <functional>
#include <set>
#include <deque>
#include <vector>
#include <string>
#include <stack>

#include <ctime>
#include <limits>
#include <iostream>

#include "trie/trie.hpp"
#include "annotation/annotation.hpp"
#include "kstat/kstat.hpp"
#include "alicont/alicont.hpp"

#include "contig_const_iterator.hpp"
#include "contig_iterator.hpp"

template <class T, template <class> class Property, class LabelType=std::string>
class contig
{
public:
    friend class contig_const_iterator<T, Property, LabelType>;
    friend class contig_iterator<T, Property, LabelType>;

    typedef Property<T>                        			  data_type;
    typedef std::vector<byte>				   			  alphabet_type;
    typedef annotation<T, Property, LabelType> 			  anno_type;
    typedef typename anno_type::index          			  link_type;
    typedef typename anno_type::record_type    			  record_type;
    typedef trie<std::vector<link_type>>				  trie_type;
    typedef typename trie_type::index_type				  index_type;
    typedef kstatistics<index_type>		      			  kstat_type;

    typedef typename trie_type::const_iterator      	  const_iterator;
    typedef typename trie_type::iterator       			  iterator;

    typedef contig_const_iterator<T, Property, LabelType> const_aiterator;
    typedef contig_iterator<T, Property, LabelType>       aiterator;

    typedef std::pair<record_type,
                      std::vector<record_type>>           alignment_type;

    contig(std::string const & name, alphabet_type const & alphabet, size_t k = 7)
        : m_name(name), m_stat(alphabet, k)
    {
    }

    template <class Iterator1, class Iterator2>
    iterator push(Iterator1 begin, Iterator1 end, Iterator2 anno_begin,
                  LabelType label)
    {
        for (Iterator1 iter = begin; iter != end; ++iter)
        {
            if (!m_stat.inAlphabet(*iter))
                return this->end();
        }

        size_t record_length = std::distance(begin, end);

        record_type & new_record = m_anno.add(label, record_length);
        iterator trie_iter = m_trie.begin();

        for (Iterator1 iter = begin; iter != end; ++iter, ++anno_begin)
        {
            byte symbol = *iter;
            data_type new_data = *anno_begin;

            // Add trie node
            trie_iter = m_trie.insert(trie_iter, symbol);

            // Add annotation node and link annotation to trie node
            data_type & data = new_record.push(trie_iter.index(), symbol);
            // Add annotation data
            data = new_data;

            // Link trie node and annotation data
            trie_iter->push_back(m_anno.last());

            // Add statistics data and link it with the node
            m_stat.add(iter, end, trie_iter.index());
        }

        return trie_iter - (record_length - 1);
    }

    template <class Iterator>
    iterator push(Iterator begin, Iterator end, LabelType label)
    {
        struct abi_false
        {
            data_type operator*() const
            {
                return data_type();
            }

            abi_false & operator++()
            {
                return *this;
            }

            abi_false & operator++(int)
            {
                return *this;
            }
        };

        return push(begin, end, abi_false(), label);
    }

    std::vector<data_type> getAnnotations(const_iterator iter) const
    {
        std::vector<data_type> result;
        for (auto i = iter->begin(); i != iter->end(); ++i)
        {
            result.push_back(m_anno[*i]);
        }
        return result;
    }

    std::vector<LabelType> getLabels(const_iterator iter) const
    {
        std::vector<LabelType> result;
        for (auto i = iter->begin(); i != iter->end(); ++i)
        {
            result.push_back(m_anno.labelOf(*i));
        }
        return result;
    }

    // Trie iterators
    iterator begin()
    {
        return m_trie.begin();
    }

    iterator end()
    {
        return m_trie.end();
    }

    const_iterator begin() const
    {
        return m_trie.begin();
    }

    const_iterator end() const
    {
        return m_trie.end();
    }

    iterator iter(index_type i)
    {
        return iterator(&m_trie, i);
    }

    const_iterator iter(index_type i) const
    {
        return const_iterator(&m_trie, i);
    }

    // Annotation iterators
    aiterator abegin()
    {
        return aiterator(this, begin());
    }

    const_aiterator abegin() const
    {
        return const_aiterator(this, begin());
    }

    aiterator aend()
    {
        return aiterator(this, end());
    }

    const_aiterator aend() const
    {
        return const_aiterator(this, end());
    }

    aiterator aiter(index_type i)
    {
        return aiterator(this, i);
    }

    const_aiterator aiter(index_type i) const
    {
        return const_aiterator(this, i);
    }

    template <class S>
    void copyTrie(trie<S> ** t) const
    {
        *t = new trie<S>(m_trie);
    }

    record_type & getRecord(LabelType const & label)
    {
        return m_anno.getRecordByLabel(label);
    }

    record_type & getRecord(size_t record)
    {
        return m_anno[record];
    }

    size_t size()
    {
        return m_trie.size();
    }

// API Methods
    template <class Iterator>
    std::vector<index_type> find(Iterator begin, Iterator end)
    {
        typedef typename trie<bool>::iterator tb_iterator;

        std::vector<index_type> result;
        trie<bool>* tmp = nullptr;
        copyTrie(&tmp);

        size_t deep = 0;
        const std::set<index_type>* nodes = nullptr;
        for (Iterator iter = begin; iter != end; ++iter)
        {
            nodes = m_stat.get(iter, end);
            if (nodes == nullptr)
            {
                nodes = m_stat.get(iter-1, end);
                break;
            }
            for (auto niter = nodes->begin(); niter != nodes->end(); ++niter)
            {
                tb_iterator t(tmp, *niter);
                *t = true;
            }
            deep++;
        }

        if (nodes == nullptr)
        {
            return result;
        }

        for (auto niter = nodes->begin(); niter != nodes->end(); ++niter)
        {
            bool breakflag = false;
            tb_iterator t(tmp, *niter);
            for (size_t i = deep-1; i != 0; --i)
            {
                --t;
                if (t == tmp->end() || !*t)
                {
                    breakflag = true;
                    break;
                }
            }
            if (!breakflag)
            {
                result.push_back(t.index());
            }
        }

        delete tmp;
        return result;
    }

    template <class Iterator>
    std::vector<alignment_result> align(Iterator begin, Iterator end,
                                        int gap, score_matrix const & matrix,
                                        size_t count)
    {
        struct results_set
        {
            std::set<alignment_result> m_set;
            size_t                     m_count;

            results_set(size_t c) : m_count(c)
            {
            }

            bool operator()(alignment_result && res)
            {
                m_set.insert(res);
                if (m_set.size() > m_count)
                {
                    m_set.erase(m_set.begin());
                }
                return true;
            }
        };

        results_set rs(count);

        align_template(begin, end, gap, matrix, rs);

        std::vector<alignment_result> results(rs.m_set.begin(), rs.m_set.end());
        return results;
    }

    template <class Iterator>
    std::vector<alignment_result> align_sc(Iterator begin, Iterator end,
                                           int gap, score_matrix const & matrix,
                                           double score_part)
    {
        struct results_vector
        {
            std::vector<alignment_result> m_vec;
            double                        m_score;

            results_vector(double s) : m_score(s)
            {
            }

            bool operator()(alignment_result && res)
            {
                if (res.score >= m_score)
                {
                    m_vec.push_back(res);
                }
                return true;
            }
        };

        std::string query(begin, end);
        results_vector rv(score_part *
                          needleman_wunsch(query, query, gap,
                                           matrix).back().back());

        align_template(begin, end, gap, matrix, rv);

        return rv.m_vec;
    }

    template <class Iterator>
    record_type annotate(Iterator begin, Iterator end, int gap,
                         score_matrix const & matrix)
    {
        std::vector<alignment_result> vec(std::move(align(begin, end, gap, matrix, 1)));
        std::string & query = vec[0].first;
        std::string & target = vec[0].second;
        size_t target_id = vec[0].target_id;
        record_type result(std::distance(begin, end));
        for (size_t i = 0; i < target.size(); ++i)
        {
            if (target[i] == '-' && i != 0)
            {
                data_type & data = result.push(std::numeric_limits<size_t>::max(), query[i]);
                data = result[result.size() - 2];
            }
            else if (target[i] == '-' && i == 0)
            {
                result.push(std::numeric_limits<size_t>::max(), query[i]);
            }
            else if (query[i] == '-')
            {
                continue;
            }
            else
            {
                data_type & data = result.push(std::numeric_limits<size_t>::max(), query[i]);
                data = getRecord(target_id)[i];
            }
        }

        return result;
    }

private:
    template <class Iterator, class Callable>
    void align_template(Iterator begin, Iterator end,
                        int gap, score_matrix const & matrix,
                        Callable & callfunc)
    {
        typedef std::pair<std::string, simple_matrix2i> node_cache_type;

        std::string query(begin, end);
        alicont ali(query, gap, matrix);
        trie<node_cache_type>* tmp;
        copyTrie(&tmp);

        std::string target;
        std::stack<size_t> index_stack;

        for (trie<node_cache_type>::iterator i = tmp->begin() + 1;
                                             i != tmp->end();
                                             ++i)
        {
            target.push_back(i.symbol());
            if (!index_stack.empty() && i.prev().fork() &&
                    i.prev().index() != index_stack.top())
            {
                while (!index_stack.empty() && index_stack.top() != i.prev().index())
                {
                    *(tmp->iter(index_stack.top())) = node_cache_type();
                    index_stack.pop();
                    ali.pop();
                }
            }
            if (i.leaf())
            {
                *i = std::move(std::make_pair(target, ali.score(target)));
                ali.push(i->first, &(i->second));
                index_stack.push(i.index());
                target.clear();
                alignment_result res = ali.alignment();
                res.target_id = iter(i.index())->back().record;
                if (!callfunc(std::move(res)))
                {
                    break;
                }
            }
            else if (i.fork())
            {
                *i = std::move(std::make_pair(target, ali.score(target)));
                ali.push(i->first, &(i->second));
                index_stack.push(i.index());
                target.clear();
            }
        }

        delete tmp;
    }

//    template <class Iterator>
//    std::vector<alignment_result> align(Iterator begin, Iterator end,
//                                        int gap, score_matrix const & matrix,
//                                        size_t count)
//    {
//        typedef std::pair<std::string, simple_matrix2i> node_cache_type;

//        std::string query(begin, end);
//        alicont ali(query, gap, matrix);
//        trie<node_cache_type>* tmp;
//        copyTrie(&tmp);

//        std::set<alignment_result> results_set;
//        std::string target;
//        std::stack<size_t> index_stack;

//        for (trie<node_cache_type>::iterator i = tmp->begin() + 1;
//                                             i != tmp->end();
//                                             ++i)
//        {
//            target.push_back(i.symbol());
//            if (!index_stack.empty() && i.prev().fork() &&
//                    i.prev().index() != index_stack.top())
//            {
//                while (!index_stack.empty() && index_stack.top() != i.prev().index())
//                {
//                    *(tmp->iter(index_stack.top())) = node_cache_type();
//                    index_stack.pop();
//                    ali.pop();
//                }
//            }
//            if (i.leaf())
//            {
//                *i = std::move(std::make_pair(target, ali.score(target)));
//                ali.push(i->first, &(i->second));
//                index_stack.push(i.index());
//                target.clear();
//                alignment_result res = ali.alignment();
//                res.target_id = iter(i.index())->back().record;
//                results_set.insert(std::move(res));
//                if (results_set.size() > count)
//                {
//                    results_set.erase(results_set.begin());
//                }
//            }
//            else if (i.fork())
//            {
//                *i = std::move(std::make_pair(target, ali.score(target)));
//                ali.push(i->first, &(i->second));
//                index_stack.push(i.index());
//                target.clear();
//            }
//        }

//        std::vector<alignment_result> results(results_set.begin(), results_set.end());

//        delete tmp;
//        return results;
//    }

private:
    std::string m_name;
    trie_type   m_trie;
    anno_type   m_anno;
    kstat_type  m_stat;
};
