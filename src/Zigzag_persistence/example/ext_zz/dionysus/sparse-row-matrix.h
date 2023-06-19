#ifndef DIONYSUS_SPARSE_ROW_MATRIX_H
#define DIONYSUS_SPARSE_ROW_MATRIX_H

#include <vector>
#include <list>
#include <unordered_map>
#include <iostream>         // for debugging output

#include <boost/intrusive/list.hpp>

#include "chain.h"
#include "reduction.h"

namespace dionysus
{

namespace bi = boost::intrusive;

namespace detail
{
    typedef         bi::list_base_hook<bi::link_mode<bi::auto_unlink>>      auto_unlink_hook;

    template<class F, class I>
    struct SparseRowMatrixEntry:
        public ChainEntry<F, std::tuple<I,I>, auto_unlink_hook>
    {
        typedef             I                                                   Index;
        typedef             typename F::Element                                 FieldElement;
        typedef             std::tuple<Index, Index>                            IndexPair;                  // (id, pair)
        typedef             ChainEntry<F, IndexPair, auto_unlink_hook>          Parent;
        typedef             SparseRowMatrixEntry                                Entry;

                            SparseRowMatrixEntry(FieldElement e, const IndexPair& ip):
                                Parent(e,ip)                                    {}

                            SparseRowMatrixEntry(FieldElement e, const Index& r, const Index& c):
                                Parent(e,IndexPair(r,c))                        {}

                            SparseRowMatrixEntry(const Entry& other)   = default;
                            SparseRowMatrixEntry(Entry&& other)        = default;
        Entry&              operator=(Entry&& other)    = default;

        void                unlink()                                            { auto_unlink_hook::unlink(); }
        bool                is_linked()  const                                  { return auto_unlink_hook::is_linked();  }
    };
}

template<class Field_, class Index_ = int, class Comparison_ = std::less<Index_>,
         template<class E, class... Args> class Column_ = std::vector>
class SparseRowMatrix
{
    public:
        typedef         Field_                                                  Field;
        typedef         Index_                                                  Index;
        typedef         Comparison_                                             Comparison;

        typedef         typename Field::Element                                 FieldElement;

        typedef         detail::SparseRowMatrixEntry<Field,Index>               Entry;
        typedef         Column_<Entry>                                          Column;
        typedef         typename Entry::IndexPair                               IndexPair;
        typedef         bi::list<Entry, bi::constant_time_size<false>>          Row;

        typedef         std::vector<ChainEntry<Field, Index>>                   IndexChain;

        typedef         std::unordered_map<Index, Column>                       Columns;
        typedef         std::unordered_map<Index, Row>                          Rows;
        typedef         std::unordered_map<Index, Index>                        LowMap;

    public:
                        SparseRowMatrix(const Field&        field,
                                        const Comparison&   cmp = Comparison()):
                            field_(field), cmp_(cmp)                            {}

                        SparseRowMatrix(SparseRowMatrix&& other)                = default;


        template<class ChainRange>
        Column          reduce(const ChainRange& chain, IndexChain& trail);

        Index           set(Index i, Column&& chain);                           // returns previous column with this low
        void            fix(Index c, Column& column);
        void            fix(Index c)                                            { fix(c, col(c)); }

        const Row&      prepend_row(Index r, FieldElement m, const Row& chain); // could be horribly inefficient if Column is chosen poorly

        void            drop_row(Index r)                                       { rows_.erase(r); if (is_low(r)) lows_.erase(r); }
        void            drop_col(Index c)
        {
            auto cit = columns_.find(c);
            Column& column = cit->second;
            if (!column.empty())
            {
                Index rlow = std::get<0>(column.back().index());
                auto it = lows_.find(rlow);
                if (it != lows_.end() && it->second == c)
                    lows_.erase(it);
            }
            columns_.erase(cit);
        }
        void            drop_low(Index r)                                       { lows_.erase(r); }

        // accessors
        Row&            row(Index r)                                            { return rows_[r]; }
        Column&         col(Index c)                                            { assert(col_exists(c)); return columns_.find(c)->second; }
        const Column&   col(Index c) const                                      { assert(col_exists(c)); return columns_.find(c)->second; }
        Index           low(Index r) const                                      { return lows_.find(r)->second; }
        bool            is_low(Index r) const                                   { return lows_.find(r) != lows_.end(); }
        void            update_low(Index c)                                     { lows_[std::get<0>(col(c).back().index())] = c; }

        const Field&        field() const                                       { return field_; }
        void                reserve(size_t)                                     {}                              // here for compatibility only
        const Comparison&   cmp() const                                         { return cmp_; }

        // debug
        bool            col_exists(Index c) const                               { return columns_.find(c) != columns_.end(); }
        const Columns&  columns() const                                         { return columns_; }
        void            check_columns() const
        {
            for (auto& x : columns_)
            {
                Index c  = x.first;
                if (x.second.empty())
                    std::cout << "Warning: empty column " << c << std::endl;
                Index rl = std::get<0>(x.second.back().index());
                if (!is_low(rl) || low(rl) != c)
                {
                    std::cout << "Columns don't check out: lows don't match" << std::endl;
                    std::cout << "   " << c << ' ' << rl << ' ' << ' ' << low(rl) << std::endl;
                    std::cout << "---\n";
                    for (auto& x : col(c))
                        std::cout << "   " << x.element() << ' ' << std::get<0>(x.index()) << ' ' << std::get<1>(x.index()) << '\n';
                    std::cout << "---\n";
                    for (auto& x : col(low(rl)))
                        std::cout << "   " << x.element() << ' ' << std::get<0>(x.index()) << ' ' << std::get<1>(x.index()) << '\n';
                    assert(0);
                }

                for (auto& x : lows_)
                {
                    if (!col_exists(x.second))
                    {
                        std::cout << "Still keeping low of a removed column" << std::endl;
                        assert(0);
                    }
                    else if (std::get<0>(col(x.second).back().index()) != x.first)
                    {
                        std::cout << "Low mismatch: " << x.second << ' ' << std::get<0>(col(x.second).back().index()) << ' ' << x.first << '\n';
                        assert(0);
                    }
                }
            }
        }

    private:
        Field           field_;
        Comparison      cmp_;

        Columns         columns_;
        Rows            rows_;
        LowMap          lows_;          // column that has this low
};


namespace detail
{

template<class Index>
struct Unpaired<std::tuple<Index,Index>>
{
    static
    constexpr std::tuple<Index,Index>
    value()
    { return std::make_tuple(std::numeric_limits<Index>::max(),
                             std::numeric_limits<Index>::max()); }
};

}

}

#include "sparse-row-matrix.hpp"

#endif
