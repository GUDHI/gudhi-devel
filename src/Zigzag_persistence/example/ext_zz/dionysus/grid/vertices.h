#ifndef VERTICES_H
#define VERTICES_H

#include <boost/iterator/iterator_facade.hpp>

namespace grid
{

template<class Vertex_>
class VerticesIterator:
    public boost::iterator_facade<VerticesIterator<Vertex_>,
                                  Vertex_,
                                  boost::forward_traversal_tag,
                                  Vertex_,
                                  std::ptrdiff_t>
{
    typedef     boost::iterator_facade<VerticesIterator,
                                       Vertex_,
                                       boost::forward_traversal_tag,
                                       Vertex_,
                                       std::ptrdiff_t>              Parent;


    public:
        typedef     typename Parent::value_type                     value_type;
        typedef     typename Parent::difference_type                difference_type;
        typedef     typename Parent::reference                      reference;

        typedef     value_type                                      Vertex;
        typedef     typename Vertex::Coordinate                     Coordinate;

                    // upper bounds are non-inclusive
                    VerticesIterator(const Vertex& bounds):
                        to_(bounds - Vertex::one())                 {}

                    VerticesIterator(const Vertex& pos,
                                     const Vertex& bounds):
                        pos_(pos), to_(bounds - Vertex::one())      {}

                    VerticesIterator(const Vertex& pos,
                                     const Vertex& from,
                                     const Vertex& to):
                        pos_(pos), from_(from),
                        to_(to)                                     {}


        static VerticesIterator
                    begin(const Vertex& bounds)                     { return VerticesIterator(bounds); }
        static VerticesIterator
                    end(const Vertex& bounds)                       { Vertex e; e[0] = bounds[0]; return VerticesIterator(e, bounds); }

        static VerticesIterator
                    begin(const Vertex& from, const Vertex& to)     { return VerticesIterator(from, from, to); }
        static VerticesIterator
                    end(const Vertex& from, const Vertex& to)       { Vertex e = from; e[0] = to[0] + 1; return VerticesIterator(e, from, to); }

    private:
        void        increment();
        bool        equal(const VerticesIterator& other) const      { return pos_ == other.pos_; }
        reference   dereference() const                             { return pos_; }

        friend class ::boost::iterator_core_access;

    private:
        Vertex      pos_;
        Vertex      from_;
        Vertex      to_;
};

}

template<class V>
void
grid::VerticesIterator<V>::
increment()
{
    unsigned j = Vertex::dimension() - 1;
    while (j > 0 && pos_[j] == to_[j])
    {
        pos_[j] = from_[j];
        --j;
    }
    ++pos_[j];
}

#endif
