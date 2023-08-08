#ifndef DIONYSUS_Z2_H
#define DIONYSUS_Z2_H

namespace dionysus
{

class Z2Field
{
    public:
        typedef         short                               Element;

                        Z2Field()                           {}

        static Element  id()                                { return 1; }
        static Element  zero()                              { return 0; }
        static Element  init(int a)                         { return (a % 2 + 2) % 2; }

        Element         neg(Element a) const                { return 2 - a; }
        Element         add(Element a, Element b) const     { return (a+b) % 2; }

        Element         inv(Element a) const                { return a; }
        Element         mul(Element a, Element b) const     { return a*b; }
        Element         div(Element a, Element b) const     { return a; }

        bool            is_zero(Element a) const            { return a == 0; }
};

}

#endif

