#ifndef DIONYSUS_CHAIN_H
#define DIONYSUS_CHAIN_H

#include <set>
#include <vector>
#include <list>

#include "fields/z2.h"

namespace dionysus
{

template<class Field>
struct FieldElement
{
    typedef     typename Field::Element     Element;
                FieldElement(Element e_):
                    e(e_)                   {}
    Element     element() const             { return e; }
    void        set_element(Element e_)     { e = e_; }
    Element     e;
};

template<>
struct FieldElement<Z2Field>
{
    typedef     Z2Field::Element            Element;
                FieldElement(Element)       {}
    Element     element() const             { return Z2Field::id(); }
    void        set_element(Element)        {}
};

template<class Field_, class Index_, class... Extra>
struct ChainEntry: public FieldElement<Field_>, public Extra...
{
    typedef     Field_                      Field;
    typedef     Index_                      Index;

    typedef     FieldElement<Field>         Parent;
    typedef     typename Parent::Element    Element;

                ChainEntry(): Parent(Element()), i(Index()) {}      // need for serialization

                ChainEntry(ChainEntry&& other)      = default;
                ChainEntry(const ChainEntry& other) = default;
    ChainEntry& operator=(ChainEntry&& other)       = default;

                ChainEntry(Element e_, const Index& i_):
                    Parent(e_), i(i_)       {}

                ChainEntry(Element e_, Index&& i_):
                    Parent(e_), i(std::move(i_))       {}

    const Index& index() const              { return i; }
    Index&      index()                     { return i; }

    // debug
    bool        operator==(const ChainEntry& other) const       { return i == other.i; }

    Index       i;
};

template<class C1>
struct Chain
{
    struct Visitor
    {
        template<class Iter>
        void first(Iter it) const               {}

        template<class Iter>
        void second(Iter it) const              {}

        template<class Iter>
        void equal_keep(Iter it) const          {}

        template<class Iter>
        void equal_drop(Iter it) const          {}
    };

    // x += a*y
    template<class C2, class Field, class Cmp, class Visitor_ = Visitor>
    static void addto(C1& x, typename Field::Element a, const C2& y, const Field& field, const Cmp& cmp, const Visitor_& = Visitor_());
};

template<class T>
struct Chain<std::list<T>>
{
    struct Visitor
    {
        template<class Iter>
        void first(Iter it) const               {}

        template<class Iter>
        void second(Iter it) const              {}

        template<class Iter>
        void equal_keep(Iter it) const          {}

        template<class Iter>
        void equal_drop(Iter it) const          {}
    };

    // x += a*y
    template<class C2, class Field, class Cmp, class Visitor_ = Visitor>
    static void addto(std::list<T>& x, typename Field::Element a, const C2& y,
                      const Field& field, const Cmp& cmp, const Visitor_& visitor = Visitor_());
};


template<class T, class TCmp>
struct Chain<std::set<T,TCmp>>
{
    struct Visitor
    {
        template<class Iter>
        void first(Iter it) const               {}

        template<class Iter>
        void second(Iter it) const              {}

        template<class Iter>
        void equal_keep(Iter it) const          {}

        template<class Iter>
        void equal_drop(Iter it) const          {}
    };

    // x += a*y
    template<class C2, class Field, class Cmp, class Visitor_ = Visitor>
    static void addto(std::set<T,TCmp>& x, typename Field::Element a, const C2& y,
                      const Field& field, const Cmp& cmp, const Visitor_& = Visitor_());

    template<class Field, class Cmp, class Visitor_ = Visitor>
    static void addto(std::set<T,TCmp>& x, typename Field::Element a, T&& y,
                      const Field& field, const Cmp& cmp, const Visitor_& = Visitor_());
};

}

//namespace std
//{
//    template<class F, class I, class... E>
//    void swap(::dionysus::ChainEntry<F,I,E...>& x, ::dionysus::ChainEntry<F,I,E...>& y)
//    {
//        std::swap(x.e, y.e);
//        std::swap(x.i, y.i);
//    }
//}

#include "chain.hpp"

#endif
