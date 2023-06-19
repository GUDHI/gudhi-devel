#include <utility>
#include <iostream>
#include <algorithm>
#include <functional>

#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/set_algorithm.hpp>

template<class D, class S>
template<class Functor, class Iterator>
void
dionysus::Rips<D,S>::
generate(Dimension k, DistanceType max, const Functor& f, Iterator bg, Iterator end) const
{
    auto neighbor = [this, max](Vertex u, Vertex v) { return this->distances()(u,v) <= max; };

    // current      = empty
    // candidates   = everything
    VertexContainer current;
    VertexContainer candidates(bg, end);
    bron_kerbosch(current, candidates, std::prev(candidates.begin()), k, neighbor, f);
}

template<class D, class S>
template<class Functor, class Iterator>
void
dionysus::Rips<D,S>::
vertex_cofaces(IndexType v, Dimension k, DistanceType max, const Functor& f, Iterator bg, Iterator end) const
{
    auto neighbor = [this, max](Vertex u, Vertex v) { return this->distances()(u,v) <= max; };

    // current      = [v]
    // candidates   = everything - [v]
    VertexContainer current; current.push_back(v);
    VertexContainer candidates;
    for (Iterator cur = bg; cur != end; ++cur)
        if (*cur != v && neighbor(v, *cur))
            candidates.push_back(*cur);

    bron_kerbosch(current, candidates, std::prev(candidates.begin()), k, neighbor, f);
}

template<class D, class S>
template<class Functor, class Iterator>
void
dionysus::Rips<D,S>::
edge_cofaces(IndexType u, IndexType v, Dimension k, DistanceType max, const Functor& f, Iterator bg, Iterator end) const
{
    auto neighbor = [this, max](Vertex u, Vertex v) { return this->distances()(u,v) <= max; };

    // current      = [u,v]
    // candidates   = everything - [u,v]
    VertexContainer current; current.push_back(u); current.push_back(v);

    VertexContainer candidates;
    for (Iterator cur = bg; cur != end; ++cur)
        if (*cur != u && *cur != v && neighbor(v,*cur) && neighbor(u,*cur))
            candidates.push_back(*cur);

    bron_kerbosch(current, candidates, std::prev(candidates.begin()), k, neighbor, f);
}

template<class D, class S>
template<class Functor, class Iterator>
void
dionysus::Rips<D,S>::
cofaces(const Simplex& s, Dimension k, DistanceType max, const Functor& f, Iterator bg, Iterator end) const
{
    namespace ba = boost::adaptors;

    auto neighbor = [this, max](Vertex u, Vertex v) { return this->distances()(u,v) <= max; };

    // current      = s
    VertexContainer current(s.begin(), s.end());

    // candidates   = everything - s     that is a neighbor of every vertex in the simplex
    VertexContainer candidates;
    boost::set_difference(std::make_pair(bg, end) |
                                ba::filtered([this,&s,&neighbor](Vertex cur)
                                             { for (auto& v : s)
                                                   if (!neighbor(v, cur))
                                                       return false;
                                             }),
                          s,
                          std::back_inserter(candidates));

    bron_kerbosch(current, candidates, std::prev(candidates.begin()), k, neighbor, f, false);
}


template<class D, class S>
template<class Functor, class NeighborTest>
void
dionysus::Rips<D,S>::
bron_kerbosch(VertexContainer&                          current,
              const VertexContainer&                    candidates,
              typename VertexContainer::const_iterator  excluded,
              Dimension                                 max_dim,
              const NeighborTest&                       neighbor,
              const Functor&                            functor,
              bool                                      check_initial)
{
    if (check_initial && !current.empty())
        functor(Simplex(current));

    if (current.size() == static_cast<size_t>(max_dim) + 1)
        return;

    for (auto cur = std::next(excluded); cur != candidates.end(); ++cur)
    {
        current.push_back(*cur);

        VertexContainer new_candidates;
        for (auto ccur = candidates.begin(); ccur != cur; ++ccur)
            if (neighbor(*ccur, *cur))
                new_candidates.push_back(*ccur);
        size_t ex = new_candidates.size();
        for (auto ccur = std::next(cur); ccur != candidates.end(); ++ccur)
            if (neighbor(*ccur, *cur))
                new_candidates.push_back(*ccur);
        excluded  = new_candidates.begin() + (ex - 1);

        bron_kerbosch(current, new_candidates, excluded, max_dim, neighbor, functor);
        current.pop_back();
    }
}

template<class Distances_, class Simplex_>
typename dionysus::Rips<Distances_, Simplex_>::DistanceType
dionysus::Rips<Distances_, Simplex_>::
distance(const Simplex& s1, const Simplex& s2) const
{
    DistanceType mx = 0;
    for (auto a : s1)
        for (auto b : s2)
            mx = std::max(mx, distances_(a,b));
    return mx;
}

template<class Distances_, class Simplex_>
typename dionysus::Rips<Distances_, Simplex_>::DistanceType
dionysus::Rips<Distances_, Simplex_>::
max_distance() const
{
    DistanceType mx = 0;
    for (IndexType a = distances_.begin(); a != distances_.end(); ++a)
        for (IndexType b = std::next(a); b != distances_.end(); ++b)
            mx = std::max(mx, distances_(a,b));
    return mx;
}

template<class Distances_, class Simplex_>
typename dionysus::Rips<Distances_, Simplex_>::DistanceType
dionysus::Rips<Distances_, Simplex_>::Evaluator::
operator()(const Simplex& s) const
{
    DistanceType mx = 0;
    for (auto a = s.begin(); a != s.end(); ++a)
        for (auto b = std::next(a); b != s.end(); ++b)
            mx = std::max(mx, distances_(*a,*b));
    return mx;
}
