#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include <boost/array.hpp>
#include <boost/operators.hpp>
#include <boost/foreach.hpp>

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

namespace grid
{

template<class Coordinate_, unsigned D>
class Point: public boost::array<Coordinate_, D>,
             private   boost::addable<          Point<Coordinate_,D>                // Point + Point
                     , boost::subtractable<     Point<Coordinate_,D>                // Point - Point
                     , boost::dividable2<       Point<Coordinate_,D>, Coordinate_   // Point / Coordinate
                     , boost::multipliable2<    Point<Coordinate_,D>, Coordinate_   // Point * Coordinate, Coordinate * Point
                     > > > >
{
    public:
        typedef             Coordinate_                             Coordinate;
        typedef             boost::array<Coordinate, D>             ArrayParent;

        typedef             Point<Coordinate, D-1>                  LPoint;
        typedef             Point<Coordinate, D+1>                  UPoint;

        template<class U>
        struct rebind       { typedef Point<U,D> type; };

    public:
                            Point()                                 { for (unsigned i = 0; i < D; ++i) (*this)[i] = 0; }
                            Point(const ArrayParent& a):
                                ArrayParent(a)                      {}
        template<class T>   Point(const Point<T, D>& p)             { for (size_t i = 0; i < D; ++i) (*this)[i] = p[i]; }
        template<class T>   Point(const T* a)                       { for (unsigned i = 0; i < D; ++i) (*this)[i] = a[i]; }
        template<class T>   Point(const std::vector<T>& a)          { for (unsigned i = 0; i < D; ++i) (*this)[i] = a[i]; }

        static
        unsigned            dimension()                             { return D; }

        static Point        zero()                                  { return Point(); }
        static Point        one()                                   { Point p; for (unsigned i = 0; i < D; ++i) p[i] = 1; return p; }

        LPoint              drop(int dim) const                     { LPoint p; unsigned c = 0; for (unsigned i = 0; i < D;   ++i) { if (i == dim) continue; p[c++] = (*this)[i]; } return p; }
        UPoint              lift(int dim, Coordinate x) const       { UPoint p; for (unsigned i = 0; i < D+1; ++i) { if (i < dim) p[i] = (*this)[i]; else if (i == dim) p[i] = x; else if (i > dim) p[i] = (*this)[i-1]; } return p; }

        using ArrayParent::operator[];

        Point&              operator+=(const Point& y)              { for (unsigned i = 0; i < D; ++i) (*this)[i] += y[i];  return *this; }
        Point&              operator-=(const Point& y)              { for (unsigned i = 0; i < D; ++i) (*this)[i] -= y[i];  return *this; }
        Point&              operator*=(Coordinate a)                { for (unsigned i = 0; i < D; ++i) (*this)[i] *= a;     return *this; }
        Point&              operator/=(Coordinate a)                { for (unsigned i = 0; i < D; ++i) (*this)[i] /= a;     return *this; }

        Point               operator-() const                       { Point res; for (unsigned i = 0; i < D; ++i) res[i] = -(*this)[i]; return res; }

        Coordinate          norm() const                            { return (*this)*(*this); }

        std::ostream&       operator<<(std::ostream& out) const     { out << (*this)[0]; for (unsigned i = 1; i < D; ++i) out << " " << (*this)[i]; return out; }
        std::istream&       operator>>(std::istream& in);

        friend
        Coordinate          operator*(const Point& x, const Point& y)   { Coordinate n = 0; for (size_t i = 0; i < D; ++i) n += x[i] * y[i]; return n; }

        template<class T>
        friend
        Coordinate          operator*(const Point<T,D>& x, const Point& y)   { Coordinate n = 0; for (size_t i = 0; i < D; ++i) n += x[i] * y[i]; return n; }

    private:
            friend class boost::serialization::access;

            template<class Archive>
            void            serialize(Archive& ar, const unsigned int version)      { ar & boost::serialization::base_object<ArrayParent>(*this); }
};

template<class C, unsigned D>
std::istream&
Point<C,D>::
operator>>(std::istream& in)
{
    std::string point_str;
    in >> point_str;        // read until ' '
    std::stringstream ps(point_str);

    char x;
    for (unsigned i = 0; i < dimension(); ++i)
    {
        ps >> (*this)[i];
        ps >> x;
    }

    return in;
}


template<class Coordinate, unsigned D>
Coordinate norm2(const Point<Coordinate,D>& p)
{ Coordinate res = 0; for (unsigned i = 0; i < D; ++i) res += p[i]*p[i]; return res; }

template<class C, unsigned D>
std::ostream&
operator<<(std::ostream& out, const Point<C,D>& p)
{ return p.operator<<(out); }

template<class C, unsigned D>
std::istream&
operator>>(std::istream& in, Point<C,D>& p)
{ return p.operator>>(in); }

}

namespace opts
{
    template<class T>
    struct Traits;

    template<class C, unsigned D>
    struct Traits< grid::Point<C,D> >
    {
        static
        std::string     type_string()           { return "POINT"; }
    };
}


#endif // POINT_H
