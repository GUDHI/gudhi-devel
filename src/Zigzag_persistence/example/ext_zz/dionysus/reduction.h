#ifndef DIONYSUS_REDUCTION_H
#define DIONYSUS_REDUCTION_H

#include <vector>
#include <tuple>
#include <functional>
#include <limits>
#include "chain.h"

namespace dionysus
{

namespace detail
{

template<class Index>
struct Unpaired
{ static constexpr Index value()          { return std::numeric_limits<Index>::max(); } };

}

template<class Index_>
struct Reduction
{
    typedef     Index_              Index;

    template<class Field>
    using AddtoVisitor = std::function<void(typename Field::Element, Index)>;

    template<class Item>
    struct CallToSub;

    static const Index unpaired;

    template<class Chain1,
             class ChainsLookup,
             class LowLookup,
             class Field,
             class Comparison = std::less<Index>>
    static
    Index reduce(Chain1&                     c,
                 const ChainsLookup&         chains,
                 const LowLookup&            lows,
                 const Field&                field,
                 const AddtoVisitor<Field>&  visitor = [](typename Field::Element, Index) {},
                 const Comparison&           cmp     = Comparison())
    {
        typedef     typename Field::Element         FieldElement;

        while (!c.empty())
        {
            //auto&  low = c.back();
            auto&  low = *(std::prev(c.end()));
            Index  l   = low.index();
            Index  cl  = lows(l);
			// std::cout << "idx: " << std::get<0>(cl) << ", " << std::get<1>(cl) << "\n";
            if (cl == unpaired)
                return l;
            else
            {
                // Reduce further
                auto&           co     = chains(cl);
                auto&           co_low = co.back();
                FieldElement    m      = field.neg(field.div(low.element(), co_low.element()));
                // c += m*co
                Chain<Chain1>::addto(c, m, co, field, cmp);
                visitor(m, cl);
            }
        }
        return unpaired;
    }

    template<class Chain1,
             class Chain2,
             class Field,
             class Comparison = std::less<Index>>
    static
    Index reduce(Chain1&                     c,
                 const std::vector<Chain2>&  chains,
                 const std::vector<Index>&   lows,
                 const Field&                field,
                 const AddtoVisitor<Field>&  visitor = [](typename Field::Element, Index) {},
                 const Comparison&           cmp     = Comparison())
    {
        return reduce(c,
                      CallToSub<Chain2>(chains),
                      CallToSub<Index>(lows),
                      field, visitor, cmp);
    }

    // This is a work-around a bug in GCC (should really be a lambda function)
    template<class Item>
    struct CallToSub
    {
                                        CallToSub(const std::vector<Item>& items_):
                                            items(items_)                       {}
        const Item&                     operator()(Index i) const               { return items[i]; }
        const std::vector<Item>&        items;
    };
};


template<class Index>
const Index
Reduction<Index>::unpaired = detail::Unpaired<Index>::value();

}

#endif
