#ifndef DIONYSUS_EXAMPLES_COMMON_H
#define DIONYSUS_EXAMPLES_COMMON_H

#include <sstream>
#include <fstream>

template<class PointContainer>
void    read_points(const std::string& infilename, PointContainer& points)
{
    typedef     typename PointContainer::value_type     Point;

    std::ifstream in(infilename.c_str());
    std::string   line;
    while(std::getline(in, line))
    {
        if (line[0] == '#') continue;               // comment line in the file
        std::stringstream linestream(line);
        double x;
        points.push_back(Point());
        while (linestream >> x)
            points.back().push_back(x);
    }
}

#endif // DIONYSUS_EXAMPLES_COMMON_H
