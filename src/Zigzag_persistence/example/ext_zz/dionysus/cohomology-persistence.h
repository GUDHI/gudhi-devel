#ifndef DIONYSUS_COHOMOLOGY_PERSISTENCE_H
#define DIONYSUS_COHOMOLOGY_PERSISTENCE_H

#include <vector>
#include <list>

#include <boost/intrusive/list.hpp>
namespace bi = boost::intrusive;

#include "reduction.h"
#include "chain.h"

namespace dionysus
{

template<class Field_, class Index_ = unsigned, class Comparison_ = std::less<Index_>>
class CohomologyPersistence
{
    public:
        typedef     Field_                                                  Field;
        typedef     Index_                                                  Index;
        typedef     Comparison_                                             Comparison;

        typedef     typename Field::Element                                 FieldElement;

        typedef     bi::list_base_hook<bi::link_mode<bi::auto_unlink>>      auto_unlink_hook;
        struct      Entry;
        struct      ColumnHead;

        typedef     std::vector<Entry>                                      Column;
        typedef     bi::list<Entry, bi::constant_time_size<false>>          Row;
        typedef     std::list<ColumnHead>                                   Columns;
        typedef     typename Columns::iterator                              ColumnsIterator;
        typedef     Column                                                  Chain;

        using       IndexColumn = std::tuple<Index, Column>;

                                CohomologyPersistence(const Field& field,
                                                      const Comparison& cmp = Comparison()):
                                    field_(field), cmp_(cmp)                {}

                                CohomologyPersistence(Field&& field,
                                                      const Comparison& cmp = Comparison()):
                                    field_(std::move(field)),
                                    cmp_(cmp)                               {}

                                CohomologyPersistence(CohomologyPersistence&& other):
                                    field_(std::move(other.field_)),
                                    cmp_(std::move(other.cmp_)),
                                    columns_(std::move(other.columns_)),
                                    rows_(std::move(other.rows_))           {}

        template<class ChainRange>
        Index                   add(const ChainRange& chain);

        template<class ChainRange>
        IndexColumn             add(const ChainRange& chain, bool keep_cocycle);

        // TODO: no skip support for now
        bool                    skip(Index) const                   { return false; }
        void                    add_skip()                          {}
        void                    set_skip(Index, bool flag = true)   {}

        const Field&            field() const                               { return field_; }
        const Columns&          columns() const                             { return columns_; }
        void                    reserve(size_t s)                           { rows_.reserve(s); }

        struct AddtoVisitor;

        static const Index      unpaired()                                  { return Reduction<Index>::unpaired; }

    private:
        Field                   field_;
        Comparison              cmp_;
        Columns                 columns_;
        std::vector<Row>        rows_;
};


template<class Field, class Index, class Cmp>
struct CohomologyPersistence<Field, Index, Cmp>::ColumnHead
{
                ColumnHead(Index i): index_(i)      {}

    Index       index() const                       { return index_; }

    Index       index_;
    Column      chain;
};

template<class Field, class Index, class Cmp>
struct CohomologyPersistence<Field, Index, Cmp>::Entry:
    public ChainEntry<Field, Index, auto_unlink_hook>
{
    typedef             ChainEntry<Field, Index, auto_unlink_hook>          Parent;

                        Entry(FieldElement e, const Index& i):              // slightly dangerous
                            Parent(e,i)                                     {}

                        Entry(FieldElement e, const Index& i, ColumnsIterator it):
                            Parent(e,i), column(it)                         {}

                        Entry(const Entry& other) = default;
                        Entry(Entry&& other) = default;

    void                unlink()                                            { auto_unlink_hook::unlink(); }
    bool                is_linked()  const                                  {  return auto_unlink_hook::is_linked();  }

    ColumnsIterator     column;     // TODO: I really don't like this overhead
};

}

#include "cohomology-persistence.hpp"

#endif
