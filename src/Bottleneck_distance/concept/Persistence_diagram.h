namespace Gudhi {
namespace Bottleneck_distance {
namespace Concept {

struct Diagram_point{
    double first;
    double second;
};

struct Persistence_Diagram
{
    const_iterator<Diagram_point> cbegin() const;
    const_iterator<Diagram_point> cend() const;
};

} //namespace Concept
} //namespace Bottleneck_distance
} //namespace Gudhi
