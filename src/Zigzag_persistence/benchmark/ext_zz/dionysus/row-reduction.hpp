#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

template<class F, typename I, class C, template<class Self> class... V>
template<class Filtration, class ReportPair>
void
dionysus::RowReduction<F,I,C,V...>::
operator()(const Filtration& filtration, const ReportPair& report_pair)
{
    using Cell = typename Filtration::Cell;
    (*this)(filtration, [](const Cell&) { return false; }, report_pair, &no_progress);
}

template<class F, typename I, class C, template<class Self> class... V>
template<class Filtration, class Relative, class ReportPair, class Progress>
void
dionysus::RowReduction<F,I,C,V...>::
operator()(const Filtration& filtration, const Relative& relative, const ReportPair& report_pair, const Progress& progress)
{
    persistence_.resize(filtration.size());

    typedef     typename Persistence::Index                     Index;
    typedef     typename Persistence::FieldElement              Element;
    typedef     typename Persistence::Chain                     Chain;
    typedef     typename Filtration::Cell                       Cell;
    typedef     ChainEntry<Field, Cell>                         CellChainEntry;
    typedef     ChainEntry<Field, Index>                        ChainEntry;

    std::vector<Chain>      rows(persistence_.size());

    auto& field = persistence_.field();

    // fill the matrix
    Index i = 0;
    for(auto& c : filtration)
    {
        progress();

        if (relative(c))
        {
            persistence_.set_skip(i);
            ++i;
            continue;
        }

        persistence_.set(i, c.boundary(field) |
                                       ba::filtered([relative](const CellChainEntry& e) { return !relative(e.index()); }) |
                                       ba::transformed([this,&filtration](const CellChainEntry& e)
                                       { return ChainEntry(e.element(), filtration.index(e.index())); }));
        if (!persistence_[i].empty())
        {
            auto& x = persistence_[i].back();
            rows[x.index()].emplace_back(x.element(),i);
        }
        ++i;
    }

    auto entry_cmp = [this](const ChainEntry& e1, const ChainEntry& e2) { return this->persistence_.cmp()(e1.index(), e2.index()); };

    // reduce the matrix from the bottom up
    for (auto it = rows.rbegin(); it != rows.rend(); ++it)
    {
        auto& row = *it;
        Index r   = rows.rend() - it - 1;

        if (row.empty())
            continue;

        // add the first column to every other column
        Index    c     = row.front().index();
        Element  e     = row.front().element();
        Chain&   first = persistence_.column(c);
        for (size_t i = 1; i < row.size(); ++i)
        {
            Index       cur_idx  = row[i].index();
            Element     cur_elem = row[i].element();
            Chain&      cur      = persistence_.column(cur_idx);
            if (cur.empty())        // zeroed out by the clearing optimization
                continue;

            Element m = field.neg(field.div(cur_elem, e));
            // cur += m*first
            ::dionysus::Chain<Chain>::addto(cur, m, first, field, entry_cmp);

            // update row
            if (!cur.empty())
            {
                ChainEntry ce = cur.back();
                auto& new_row = rows[ce.index()];
                new_row.emplace_back(ce.element(), cur_idx);
                if (entry_cmp(new_row.back(), new_row.front()))
                    std::swap(new_row.back(), new_row.front());
            }
        }

        persistence_.set_pair(r,c);
        report_pair(filtration[r].dimension(), r, c);

        // zero out the corresponding column (the clearing optimization)
        persistence_.column(r).clear();
    }
}

