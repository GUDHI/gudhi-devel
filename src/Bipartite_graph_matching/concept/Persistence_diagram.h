typedef typename std::pair<double,double> Diagram_point;

struct Persistence_Diagram
{
    const_iterator<Diagram_point> cbegin() const;
    const_iterator<Diagram_point> cend() const;
};
