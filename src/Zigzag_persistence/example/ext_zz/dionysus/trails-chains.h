#ifndef DIONYSUS_TRAILS_CHAINS_H
#define DIONYSUS_TRAILS_CHAINS_H

#include "ordinary-persistence.h"

template<class Field, class Index>
struct ChainsVisitor: public EmptyVisitor<Field, Index>
{
    template<class Chain>
    void        chain_initialized(Chain& c)         {  }

    void        addto(typename Field::Element m, Index cl)          {}
    void        reduction_finished()                                {}
};


#endif
