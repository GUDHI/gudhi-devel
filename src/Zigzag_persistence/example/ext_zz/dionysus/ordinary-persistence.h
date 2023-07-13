#ifndef DIONYSUS_ORDINARY_PERSISTENCE_H
#define DIONYSUS_ORDINARY_PERSISTENCE_H

#include "reduced-matrix.h"

namespace dionysus
{

/* Move this into a ReducedMatrix class */

// Ordinary D -> R reduction
template<class    Field,
         typename Index = unsigned,
         class    Comparison = std::less<Index>,
         template<class Self> class... Visitors>
using OrdinaryPersistence = ReducedMatrix<Field, Index, Comparison, Visitors...>;

// No negative optimization
template<class Field, typename Index = unsigned, class Comparison = std::less<Index>>
struct NoNegative
{
    template<class Self>
    struct Visitor: public EmptyVisitor<Field, Index, Self>
    {
        template<class Chain>
        void        chain_initialized(Self* matrix, Chain& c)
        {
            for (auto cur = std::begin(c); cur != std::end(c); ++cur)
            {
                Index i = cur->index();
                Index p = matrix->pair(i);
                if (!(p == Self::unpaired() || (*matrix)[i].empty()))
                    c.erase(cur--);
            }
        }
    };

    template<class Self>
    using V2 = EmptyVisitor<Field, Index, Self>;
};

template<class    Field,
         typename Index = unsigned,
         class    Comparison = std::less<Index>,
         template<class Self> class... Visitors>
using OrdinaryPersistenceNoNegative = ReducedMatrix<Field, Index, Comparison,
                                                    NoNegative<Field, Index, Comparison>::template Visitor,
                                                    Visitors...>;

// TODO: add clearing optimization (possibly bake it into the code itself)

template<class    Field,
         typename Index = unsigned,
         class    Comparison = std::less<Index>,
         template<class Self> class... Visitors>
using FastPersistence = ReducedMatrix<Field, Index, Comparison,
                                      NoNegative<Field, Index, Comparison>::template Visitor,
                                      //Clearing<Field, Index, Comparison>::template Visitor,                    // FIXME
                                      Visitors...>;


}

#endif
