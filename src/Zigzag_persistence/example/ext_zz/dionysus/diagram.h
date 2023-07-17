#ifndef DIONYSUS_DIAGRAM_H
#define DIONYSUS_DIAGRAM_H

#include <vector>
#include <tuple>
#include <limits>

namespace dionysus
{

template<class Value_, class Data_>
class Diagram
{
    public:
        using Value = Value_;
        using Data  = Data_;
        struct Point: public std::pair<Value, Value>
        {
            using Parent = std::pair<Value, Value>;

                  Point(Value b, Value d, Data dd):
                      Parent(b,d), data(dd)     {}

            Value   birth() const               { return Parent::first; }
            Value   death() const               { return Parent::second; }

            // FIXME: temporary hack
            Value   operator[](size_t i) const  { if (i == 0) return birth(); return death(); }

            Data    data;
        };

        using Points         = std::vector<Point>;
        using iterator       = typename Points::iterator;
        using const_iterator = typename Points::const_iterator;
        using value_type     = Point;

    public:
        const_iterator  begin() const           { return points.begin(); }
        const_iterator  end() const             { return points.end(); }
        iterator        begin()                 { return points.begin(); }
        iterator        end()                   { return points.end(); }

        const Point&    operator[](size_t i) const  { return points[i]; }

        size_t  size() const                    { return points.size(); }
        void    push_back(const Point& p)       { points.push_back(p); }
        template<class... Args>
        void    emplace_back(Args&&... args)    { points.emplace_back(std::forward<Args>(args)...); }

    private:
        std::vector<Point>      points;
};

namespace detail
{
    template<class ReducedMatrix, class Filtration, class GetValue, class GetData>
    struct Diagrams
    {
        using Value = decltype(std::declval<GetValue>()(std::declval<typename Filtration::Cell>()));
        using Data  = decltype(std::declval<GetData>()(std::declval<typename ReducedMatrix::Index>()));
        using type  = std::vector<Diagram<Value, Data>>;
    };
}

template<class ReducedMatrix, class Filtration, class GetValue, class GetData>
typename detail::Diagrams<ReducedMatrix, Filtration, GetValue, GetData>::type
init_diagrams(const ReducedMatrix& m, const Filtration& f, const GetValue& get_value, const GetData& get_data)
{
    using Result  = typename detail::Diagrams<ReducedMatrix, Filtration, GetValue, GetData>::type;

    Result diagrams;
    for (typename ReducedMatrix::Index i = 0; i < m.size(); ++i)
    {
        if (m.skip(i))
            continue;

        auto  s = f[i];
        auto  d = s.dimension();

        while (d + 1 > diagrams.size())
            diagrams.emplace_back();

        auto pair = m.pair(i);
        if (pair == m.unpaired())
        {
            auto  birth = get_value(s);
            using Value = decltype(birth);
            Value death = std::numeric_limits<Value>::infinity();
            diagrams[d].emplace_back(birth, death, get_data(i));
        } else if (pair > i)       // positive
        {
            auto birth = get_value(s);
            auto death = get_value(f[pair]);

            // hack to work with coboundaries
            auto pd = f[pair].dimension();
            if (pd < d)
            {
                d = pd;
                std::swap(birth, death);
            }

            if (birth != death)         // skip diagonal
                diagrams[d].emplace_back(birth, death, get_data(i));
        } // else negative: do nothing
    }

    return diagrams;
}

}

#endif
