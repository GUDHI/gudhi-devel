#ifndef DIONYSUS_SIMPLEX_H
#define DIONYSUS_SIMPLEX_H

#include <algorithm>
#include <functional>

//#include <boost/compressed_pair.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/functional/hash.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include "chain.h"

namespace dionysus
{

struct Empty {};

template<class Vertex_ = unsigned, class T = Empty>
class Simplex
{
    public:
        typedef         Vertex_                                     Vertex;
        typedef         T                                           Data;
        typedef         std::unique_ptr<Vertex[]>                   Vertices;

        template<class Field>
        struct BoundaryChainIterator;
        struct BoundaryIterator;

        template<class Field>
        using BoundaryChainRange = boost::iterator_range<BoundaryChainIterator<Field>>;
        using BoundaryRange      = boost::iterator_range<BoundaryIterator>;

        template<class Field>
        using Entry = ChainEntry<Field, Simplex>;

    public:
                        Simplex(const Data& d = Data()):
                            dim_(-1), data_(d)                                                  {}

                        Simplex(const std::initializer_list<Vertex>& vertices,
                                Data&& d = Data()):
                            Simplex(vertices.size() - 1, vertices.begin(), vertices.end(), std::move(d))
                                                                                                {}

                        Simplex(const std::initializer_list<Vertex>& vertices,
                                const Data& d):
                            Simplex(vertices.size() - 1, vertices.begin(), vertices.end(), d)   {}

                        Simplex(short unsigned dim, Vertices&& vertices, Data&& data = Data()):
                            dim_(dim), vertices_(std::move(vertices)), data_(std::move(data))   { std::sort(begin(), end()); }

        template<class VertexRange>
                        Simplex(const VertexRange& vertices,
                                Data&& d = Data()):
                            Simplex(vertices.size() - 1, vertices.begin(), vertices.end(), std::move(d))
                                                                                                {}

        template<class VertexRange>
                        Simplex(const VertexRange& vertices,
                                const Data& d):
                            Simplex(vertices.size() - 1, vertices.begin(), vertices.end(), d)   {}

                        Simplex(const Simplex& other):
                            Simplex(other.dim_, other.begin(), other.end(), other.data_)        {}
        Simplex&        operator=(const Simplex& other)             { dim_ = other.dim_; vertices_ = Vertices(new Vertex[dim_+1]); std::copy(other.begin(), other.end(), begin()); data_ = other.data_; return *this; }

                        Simplex(Simplex&& other) noexcept:
                            dim_(other.dim_),
                            vertices_(std::move(other.vertices_)),
                            data_(std::move(other.data_))           {}
        Simplex&        operator=(Simplex&& other)                  = default;

        template<class Iterator>
                        Simplex(short unsigned dim,
                                Iterator b, Iterator e,
                                Data&& d = Data()):
                            dim_(dim),
                            vertices_(new Vertex[dim_+1]),
                            data_(std::move(d))                     { std::copy(b, e, begin()); std::sort(begin(), end()); }

        template<class Iterator>
                        Simplex(short unsigned dim,
                                Iterator b, Iterator e,
                                const Data& d):
                            dim_(dim),
                            vertices_(new Vertex[dim_+1]),
                            data_(d)                                { std::copy(b, e, begin()); std::sort(begin(), end()); }

        short unsigned  dimension() const                           { return dim_; }

        BoundaryRange    boundary() const                           { return BoundaryRange(boundary_begin(), boundary_end()); }
        BoundaryIterator boundary_begin() const;
        BoundaryIterator boundary_end() const;

        template<class Field>
        BoundaryChainRange<Field>
                        boundary(const Field& field) const          { return BoundaryChainRange<Field>(boundary_begin(field), boundary_end(field)); }

        template<class Field>
        BoundaryChainIterator<Field>
                        boundary_begin(const Field& field) const;
        template<class Field>
        BoundaryChainIterator<Field>
                        boundary_end(const Field& field) const;

        const Vertex*   begin() const                               { return vertices_.get(); }
        const Vertex*   end() const                                 { return begin() + dim_ + 1; }
        size_t          size() const                                { return dim_ + 1; }

        std::pair<const Vertex*, const Vertex*>
                        range() const                               { return std::make_pair(begin(), end()); }

        Simplex         join(const Vertex& v) const                 { Vertices vertices(new Vertex[dim_+2]); std::copy(begin(), end(), vertices.get()); vertices[dim_+1] = v; return Simplex(dim_ + 1, std::move(vertices), Data(data_)); }

        bool            operator==(const Simplex& other) const      { return dim_ == other.dim_ && std::equal(begin(), end(), other.begin()); }
        bool            operator!=(const Simplex& other) const      { return !operator==(other); }
        bool            operator<(const Simplex& other) const       { return dim_ < other.dim_ || (dim_ == other.dim_ && std::lexicographical_compare(begin(), end(), other.begin(), other.end())); }
        bool            operator>(const Simplex& other) const       { return other < (*this); }

        Vertex          operator[](short unsigned i) const          { return vertices_[i]; }
        const Data&     data() const                                { return data_; }
        Data&           data()                                      { return data_; }

        friend
        std::ostream&   operator<<(std::ostream& out, const Simplex& s)
        { out << '<' << *s.begin(); for (auto it = s.begin() + 1; it != s.end(); ++it) out << ',' << *it; out << '>'; return out; }

    private:
        Vertex*         begin()                                     { return vertices_.get(); }
        Vertex*         end()                                       { return begin() + dim_ + 1; }

    private:
        short unsigned      dim_;
        //boost::compressed_pair<Vertices, Data>      vertices_data_;
        Vertices            vertices_;
        Data                data_;          // TODO: optimize
};

template<class V, class D>
size_t hash_value(const Simplex<V,D>& s)                            { return boost::hash_range(s.begin(), s.end()); }


template<class V, class D>
struct Simplex<V,D>::BoundaryIterator:
    public boost::iterator_adaptor<BoundaryIterator,    // Derived
                                   const V*,            // Base
                                   Simplex<V,D>,        // Value
                                   boost::use_default,
                                   Simplex<V,D>>        // Reference
{
    public:
        typedef     const V*                            Iterator;
        typedef     Simplex<V,D>                        Value;

        typedef     boost::iterator_adaptor<BoundaryIterator,
                                            Iterator,
                                            Value,
                                            boost::use_default,
                                            Value>                          Parent;

                    BoundaryIterator()                                      {}
        explicit    BoundaryIterator(short unsigned dim, Iterator iter, Iterator bg, Iterator end):
                        Parent(iter), dim_(dim), bg_(bg), end_(end)         {}

        Iterator    begin() const                                           { return bg_; }

    private:
        friend class    boost::iterator_core_access;
        Value    dereference() const
        {
            typedef     std::not_equal_to<V>                                NotEqualVertex;

            using std::placeholders::_1;
            return      Simplex(dim_ - 1,
                                boost::make_filter_iterator(std::bind(NotEqualVertex(), _1, *(this->base())), bg_,  end_),
                                boost::make_filter_iterator(std::bind(NotEqualVertex(), _1, *(this->base())), end_, end_));
        }

        short unsigned  dim_;
        Iterator        bg_;
        Iterator        end_;
};

template<class V, class D>
template<class F>
struct Simplex<V,D>::BoundaryChainIterator:
    public boost::iterator_adaptor<BoundaryChainIterator<F>,              // Derived
                                   BoundaryIterator,
                                   ChainEntry<F,Simplex<V,D>>,    // Value
                                   boost::use_default,
                                   ChainEntry<F,Simplex<V,D>>>    // Reference
{
    public:
        typedef     F                                                       Field;
        typedef     BoundaryIterator                                        Iterator;
        typedef     ChainEntry<F,Simplex<V,D>>                              Value;

        typedef     boost::iterator_adaptor<BoundaryChainIterator,
                                            Iterator,
                                            Value,
                                            boost::use_default,
                                            Value>                          Parent;

                    BoundaryChainIterator()                                      {}
        explicit    BoundaryChainIterator(const Field& field, Iterator iter):
                        Parent(iter), field_(&field)                         {}

    private:
        friend class    boost::iterator_core_access;
        Value    dereference() const
        {
            return      Value(((this->base().base() - this->base().begin()) % 2 == 0)? field_->id() : field_->neg(field_->id()),
                              *(this->base()));
        }

        const Field*    field_ = nullptr;
};


/* Simplex */
template<class V, class D>
typename Simplex<V,D>::BoundaryIterator
Simplex<V,D>::
boundary_begin() const
{
    if (dimension() == 0)   return boundary_end();
    return BoundaryIterator(dimension(), begin(), begin(), end());
}

template<class V, class D>
typename Simplex<V,D>::BoundaryIterator
Simplex<V,D>::
boundary_end() const
{
    return BoundaryIterator(dimension(), end(), begin(), end());
}

template<class V, class D>
template<class F>
#if defined(_MSC_VER)
typename Simplex<V,D>::BoundaryChainIterator<F>
#else
typename Simplex<V,D>::template BoundaryChainIterator<F>
#endif
Simplex<V,D>::
boundary_begin(const F& field) const
{
    if (dimension() == 0)   return boundary_end(field);
    return BoundaryChainIterator<F>(field, boundary_begin());
}

template<class V, class D>
template<class F>
#if defined(_MSC_VER)
typename Simplex<V,D>::BoundaryChainIterator<F>
#else
typename Simplex<V,D>::template BoundaryChainIterator<F>
#endif
Simplex<V,D>::
boundary_end(const F& field) const
{
    return BoundaryChainIterator<F>(field, boundary_end());
}

} // dionysus

namespace std
{

template<class V, class T>
struct hash<dionysus::Simplex<V,T>>
{
    size_t operator()(const dionysus::Simplex<V,T>& s) const            { return hash_value(s); }
};

} // std

#endif
