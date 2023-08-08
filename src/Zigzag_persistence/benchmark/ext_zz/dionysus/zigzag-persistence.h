#ifndef DIONYSUS_ZIGZAG_PERSISTENCE_H
#define DIONYSUS_ZIGZAG_PERSISTENCE_H

#include <tuple>
#include <type_traits>
#include <deque>

#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/map.hpp>

#include "sparse-row-matrix.h"

namespace dionysus
{

namespace ba = boost::adaptors;

template<class Field_, class Index_ = int, class Comparison_ = std::less<Index_>>
class ZigzagPersistence
{
    static_assert(std::is_signed<Index_>::value, "Index type used in ZigzagPersistence must be a *signed* integer");

    public:
        typedef         Field_                                      Field;
        typedef         Index_                                      Index;
        typedef         Comparison_                                 Comparison;

        typedef         SparseRowMatrix<Field, Index, Comparison>   RowMatrix;
        typedef         SparseRowMatrix<Field, Index, Comparison,
                                        std::deque>                 DequeRowMatrix;
        typedef         typename RowMatrix::IndexPair               IndexPair;
        typedef         typename RowMatrix::FieldElement            FieldElement;
        typedef         typename RowMatrix::IndexChain              IndexChain;
        typedef         typename RowMatrix::Column                  Column;
        typedef         typename RowMatrix::Row                     Row;
        typedef         typename DequeRowMatrix::Column             DequeColumn;
        typedef         typename DequeRowMatrix::Row                DequeRow;

        typedef         std::unordered_map<Index, Index>            BirthIndexMap;


                        ZigzagPersistence(const Field&      field,
                                          const Comparison& cmp = Comparison()):
                            Z(field, cmp), C(field, cmp), B(field, cmp),
                            operations(0),
                            cell_indices(0),
                            z_indicies_last(0),
                            z_indicies_first(-1),
                            b_indices(0)                            {}

        template<class ChainRange>
        Index           add(const ChainRange& chain)                // returns the id of the dying cycle (or unpaired)
        {
            Index res = add_impl(chain);
#ifdef DIONYSUS_ZIGZAG_DEBUG
            check_sorted();
            check_b_cols();
            Z.check_columns();
#endif
            return res;
        }
        Index           remove(Index cell)
        {
            Index res = remove_impl(cell);
#ifdef DIONYSUS_ZIGZAG_DEBUG
            check_sorted();
            check_b_cols();
            Z.check_columns();
#endif
            return res;
        }

        struct IsAlive
        {
                    IsAlive(const ZigzagPersistence& zz_): zz(&zz_)      {}
            bool    operator()(const std::pair<Index,Index>& x) const   { return zz->is_alive(x.first); }
            const   ZigzagPersistence*  zz;
        };

        bool                is_alive(Index x) const                 { return !B.is_low(x); }

        auto                alive_ops() const -> decltype(BirthIndexMap() | ba::filtered(IsAlive(*this)) | ba::map_values)
        { return birth_index | ba::filtered(IsAlive(*this)) | ba::map_values; }

        auto                alive_cycles() const -> decltype(BirthIndexMap() | ba::filtered(IsAlive(*this)) | ba::map_keys)
        { return birth_index | ba::filtered(IsAlive(*this)) | ba::map_keys; }

        size_t              alive_size() const                      { return Z.columns().size() - B.columns().size(); }

        void                reserve(size_t)                         {}              // here for compatibility only
        const Field&        field() const                           { return Z.field(); }
        const Comparison&   cmp() const                             { return Z.cmp(); }

        template<class Entry>
        static Index    row(const Entry& e)                         { return std::get<0>(e.index()); }
        template<class Entry>
        static Index    col(const Entry& e)                         { return std::get<1>(e.index()); }

        static
        const Index     unpaired()                                  { return Reduction<Index>::unpaired; }

        const Column&   cycle(Index i) const                        { return Z.col(i); }

        // debug
        void            check_b_cols() const;

        template<class SimplexToIndex, class IndexToSimplex>
        void            check_boundaries(const SimplexToIndex& s2i, const IndexToSimplex& i2s) const;
        template<class SimplexToIndex, class IndexToSimplex>
        void            check_cycles(const SimplexToIndex& s2i, const IndexToSimplex& i2s) const;

        Column          zb_dot(Index c) const;

        template<class SimplexToIndex, class IndexToSimplex>
        Column          dc_dot(Index c, const SimplexToIndex& s2i, const IndexToSimplex& i2s) const;

        template<class SimplexToIndex, class IndexToSimplex>
        Column          boundary(Index i, const SimplexToIndex& s2i, const IndexToSimplex& i2s) const;

        void            check_sorted() const;

    private:
        template<class ChainRange>
        Index           add_impl(const ChainRange& chain);
        Index           remove_impl(Index cell);

    private:
        RowMatrix       Z, C;
        DequeRowMatrix  B;

        BirthIndexMap   birth_index;
        Index           operations;
        Index           cell_indices;
        Index           z_indicies_last, z_indicies_first;
        Index           b_indices;
};

}

#include "zigzag-persistence.hpp"

#endif
