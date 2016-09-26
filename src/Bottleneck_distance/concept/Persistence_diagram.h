namespace Gudhi {
namespace Bottleneck_distance {
namespace Concept {

/** \brief Concept of persistence diagram point. The double first is the birth of the component and the double second is the death of the component.
 *
 * \ingroup bottleneck_distance
 */
struct Diagram_point{
    double first;
    double second;
};

/** \brief Concept of persistence diagram.
 *
 * \ingroup bottleneck_distance
 */
struct Persistence_Diagram
{
    const_iterator<Diagram_point> cbegin() const;
    const_iterator<Diagram_point> cend() const;
};

} //namespace Concept
} //namespace Bottleneck_distance
} //namespace Gudhi
