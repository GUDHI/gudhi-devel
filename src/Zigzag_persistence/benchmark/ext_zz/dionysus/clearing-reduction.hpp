#include <numeric>
#include <algorithm>

#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

template<class P>
template<class Filtration, class ReportPair>
void
dionysus::ClearingReduction<P>::
operator()(const Filtration& filtration, const ReportPair& report_pair)
{
    using Cell = typename Filtration::Cell;
    (*this)(filtration, [](const Cell&) { return false; }, report_pair, &no_progress);
}

template<class P>
template<class Filtration, class Relative, class ReportPair, class Progress>
void
dionysus::ClearingReduction<P>::
operator()(const Filtration& filtration, const Relative& relative, const ReportPair& report_pair, const Progress& progress)
{
    persistence_.resize(filtration.size());

    // sort indices by decreasing dimension
    std::vector<size_t> indices(filtration.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::stable_sort(indices.begin(), indices.end(),
                     [&filtration](size_t x, size_t y)
                     { return filtration[x].dimension() > filtration[y].dimension(); });

    typedef     typename Filtration::Cell                       Cell;
    typedef     ChainEntry<Field, Cell>                         CellChainEntry;
    typedef     ChainEntry<Field, Index>                        ChainEntry;

    for(size_t i : indices)
    {
        progress();
        const auto& c = filtration[i];

        if (relative(c))
        {
            persistence_.set_skip(i);
            continue;
        }

        if (persistence_.pair(i) != persistence_.unpaired())
            continue;

        persistence_.set(i, c.boundary(persistence_.field()) |
                                       ba::filtered([relative](const CellChainEntry& e) { return !relative(e.index()); }) |
                                       ba::transformed([this,&filtration](const CellChainEntry& e)
                                       { return ChainEntry(e.element(), filtration.index(e.index())); }));

        Index pair = persistence_.reduce(i);
        if (pair != persistence_.unpaired())
            report_pair(c.dimension(), pair, i);
    }
}

