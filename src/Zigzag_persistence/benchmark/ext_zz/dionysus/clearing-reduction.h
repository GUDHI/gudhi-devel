#ifndef DIONYSUS_CLEARING_REDUCTION_H
#define DIONYSUS_CLEARING_REDUCTION_H

namespace dionysus
{

// Mid-level interface
template<class Persistence_>
class ClearingReduction
{
    public:
        using Persistence = Persistence_;
        using Field       = typename Persistence::Field;
        using Index       = typename Persistence::Index;

    public:
                    ClearingReduction(Persistence& persistence):
                        persistence_(persistence)               {}

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
        Persistence&  persistence_;
};

}

#include "clearing-reduction.hpp"

#endif

