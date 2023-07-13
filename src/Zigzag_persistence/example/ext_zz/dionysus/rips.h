#ifndef DIONYSUS_RIPS_H
#define DIONYSUS_RIPS_H

#include <vector>
#include <string>

#include <boost/iterator/counting_iterator.hpp>

#include "simplex.h"

namespace dionysus
{

/**
 * Rips class
 *
 * Class providing basic operations to work with Rips complexes. It implements Bron-Kerbosch algorithm,
 * and provides simple wrappers for various functions.
 *
 * Distances_ is expected to define types IndexType and DistanceType as well as
 *               provide operator()(...) which given two IndexTypes should return
 *               the distance between them. There should be methods begin() and end()
 *               for iterating over IndexTypes as well as a method size().
 */
template<class Distances_, class Simplex_ = Simplex<typename Distances_::IndexType> >
class Rips
{
    public:
        typedef             Distances_                                      Distances;
        typedef             typename Distances::IndexType                   IndexType;
        typedef             typename Distances::DistanceType                DistanceType;

        typedef             Simplex_                                        Simplex;
        typedef             typename Simplex::Vertex                        Vertex;             // should be the same as IndexType
        typedef             std::vector<Vertex>                             VertexContainer;

        typedef             short unsigned                                  Dimension;

        class               Evaluator;
        class               Comparison;

    public:
                            Rips(const Distances& distances):
                                distances_(distances)                       {}

        // Calls functor f on each simplex in the k-skeleton of the Rips complex
        template<class Functor, class Iterator>
        void                generate(Dimension k, DistanceType max, const Functor& f,
                                     Iterator candidates_begin, Iterator candidates_end) const;

        // Calls functor f on all the simplices of the Rips complex that contain the given vertex v
        template<class Functor, class Iterator>
        void                vertex_cofaces(IndexType v, Dimension k, DistanceType max, const Functor& f,
                                           Iterator candidates_begin, Iterator candidates_end) const;

        // Calls functor f on all the simplices of the Rips complex that contain the given edge [u,v]
        template<class Functor, class Iterator>
        void                edge_cofaces(IndexType u, IndexType v, Dimension k, DistanceType max, const Functor& f,
                                         Iterator candidates_begin, Iterator candidates_end) const;

        // Calls functor f on all the simplices of the Rips complex that contain the given Simplex s
        // (unlike the previous methods it does not call the functor on the Simplex s itself)
        template<class Functor, class Iterator>
        void                cofaces(const Simplex& s, Dimension k, DistanceType max, const Functor& f,
                                    Iterator candidates_begin, Iterator candidates_end) const;


        /* No Iterator argument means Iterator = IndexType and the range is [distances().begin(), distances().end()) */
        template<class Functor>
        void                generate(Dimension k, DistanceType max, const Functor& f) const
        { generate(k, max, f, boost::make_counting_iterator(distances().begin()), boost::make_counting_iterator(distances().end())); }

        template<class Functor>
        void                vertex_cofaces(IndexType v, Dimension k, DistanceType max, const Functor& f) const
        { vertex_cofaces(v, k, max, f, boost::make_counting_iterator(distances().begin()), boost::make_counting_iterator(distances().end())); }

        template<class Functor>
        void                edge_cofaces(IndexType u, IndexType v, Dimension k, DistanceType max, const Functor& f) const
        { edge_cofaces(u, v, k, max, f, boost::make_counting_iterator(distances().begin()), boost::make_counting_iterator(distances().end())); }

        template<class Functor>
        void                cofaces(const Simplex& s, Dimension k, DistanceType max, const Functor& f) const
        { cofaces(s, k, max, f, boost::make_counting_iterator(distances().begin()), boost::make_counting_iterator(distances().end())); }


        const Distances&    distances() const                               { return distances_; }
        DistanceType        max_distance() const;

        DistanceType        distance(const Simplex& s1, const Simplex& s2) const;


        template<class Functor, class NeighborTest>
        static void         bron_kerbosch(VertexContainer&                          current,
                                          const VertexContainer&                    candidates,
                                          typename VertexContainer::const_iterator  excluded,
                                          Dimension                                 max_dim,
                                          const NeighborTest&                       neighbor,
                                          const Functor&                            functor,
                                          bool                                      check_initial = true);

    protected:
        const Distances&    distances_;
};

template<class Distances_, class Simplex_>
class Rips<Distances_, Simplex_>::Evaluator
{
    public:
        typedef             Simplex_                                        Simplex;

                            Evaluator(const Distances& distances):
                                distances_(distances)                       {}

        DistanceType        operator()(const Simplex& s) const;

    protected:
        const Distances&    distances_;
};

template<class Distances_, class Simplex_>
class Rips<Distances_, Simplex_>::Comparison
{
    public:
        typedef             Simplex_                                        Simplex;

                            Comparison(const Distances& distances):
                                eval_(distances)                            {}

        bool                operator()(const Simplex& s1, const Simplex& s2) const
        {
            DistanceType e1 = eval_(s1),
                         e2 = eval_(s2);
            if (e1 == e2)
                return s1.dimension() < s2.dimension();

            return e1 < e2;
        }

    protected:
        Evaluator           eval_;
};

}

#include "rips.hpp"

#endif // DIONYSUS_RIPS_H
