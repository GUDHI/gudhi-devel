#ifndef DIONYSUS_ROW_REDUCTION_H
#define DIONYSUS_ROW_REDUCTION_H

#include "reduced-matrix.h"

namespace dionysus
{

// Mid-level interface
template<class Field_, typename Index_ = unsigned, class Comparison_ = std::less<Index_>, template<class Persistence> class... Visitors>
class RowReduction
{
    public:
        typedef         Field_                                                  Field;
        typedef         Index_                                                  Index;
        typedef         Comparison_                                             Comparison;

        typedef         ReducedMatrix<Field_,Index_,Comparison_,Visitors...>    Persistence;

    public:
                        RowReduction(const Field& field):
                            persistence_(field)                         {}

                        RowReduction(const Field&                       field,
                                     const Comparison&                  cmp,
                                     const Visitors<Persistence>&...    visitors):
                            persistence_(field, cmp, visitors...)       {}

        template<class Filtration, class Relative, class ReportPair, class Progress>
        void            operator()(const Filtration& f, const Relative& relative, const ReportPair& report_pair, const Progress& progress);

        template<class Filtration, class ReportPair>
        void            operator()(const Filtration& f, const ReportPair& report_pair);

        template<class Filtration>
        void            operator()(const Filtration& f)             { return (*this)(f, &no_report_pair); }

        static void     no_report_pair(int, Index, Index)           {}
        static void     no_progress()                               {}

        const Persistence&
                        persistence() const                         { return persistence_; }
        Persistence&    persistence()                               { return persistence_; }

    private:
        Persistence     persistence_;
};

}

#include "row-reduction.hpp"

#endif

