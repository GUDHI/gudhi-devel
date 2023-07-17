#ifndef RELATIVE_HOMOLOGY_ZIGZAG_H
#define RELATIVE_HOMOLOGY_ZIGZAG_H

#include <boost/bimap.hpp>
#include <boost/range/adaptors.hpp>

#include "zigzag-persistence.h"

namespace dionysus
{

namespace ba = boost::adaptors;

template<class Field_, class Index_ = int, class Comparison_ = std::less<Index_>>
class RelativeHomologyZigzag
{
    public:
        typedef         Field_                                      Field;
        typedef         Index_                                      Index;
        typedef         Comparison_                                 Comparison;

        typedef         ZigzagPersistence<Field, Index, Comparison> ZZP;
        typedef         typename ZZP::IndexChain                    IndexChain;
        typedef         typename ZZP::FieldElement                  FieldElement;
        typedef         typename IndexChain::value_type             ChainEntry;


        typedef         Comparison                                  Cmp;

                        RelativeHomologyZigzag(const Field&      field,
                                               const Comparison& cmp = Comparison()):
                            zzp_(field, cmp)
        {
            zzp_.add( IndexChain() );       // vertex w
            ++zzp_op_;
            ++zzp_cell_;
        }

        template<class ChainRange>
        void            add_both(const ChainRange& chain);

        void            remove_both(Index cell);

        // index of the absolute cell; chain = its boundary
        template<class ChainRange>
        Index           add(Index cell, const ChainRange& chain);   // add to the relative part

        Index           remove(Index cell);                         // remove from the relative part

        const Field&    field() const                               { return zzp_.field(); }
        const Cmp&      cmp() const                                 { return zzp_.cmp(); }

        size_t          alive_size() const                          { return zzp_.alive_size() - 1; }   // -1 for the cone vertex

        static
        const Index     unpaired()                                  { return ZZP::unpaired(); }

    private:
        template<class ChainRange>
        IndexChain      relative_chain(Index cell, const ChainRange& chain) const;

        template<class ChainRange>
        IndexChain      absolute_chain(const ChainRange& chain) const;

        Index           abs_index(Index idx) const                  { return absolute_.left.find(idx)->second; }
        Index           rel_index(Index idx) const                  { return relative_.left.find(idx)->second; }
        Index           decode_pair(Index pair);

    private:
        ZZP                                 zzp_;                       // underlying (cone) implementation
        boost::bimap<Index, Index>          absolute_;                  // bimap between our cells and zzp absolute cells
        boost::bimap<Index, Index>          relative_;                  // bimap between our cells and zzp relative cells
        std::unordered_map<Index, Index>    op_map_;                    // map from zzp_op to our op
        Index                               op_         = 0,
                                            zzp_op_     = 0,
                                            cell_       = 0,
                                            zzp_cell_   = 0;
};

}

#include "relative-homology-zigzag.hpp"

#endif
