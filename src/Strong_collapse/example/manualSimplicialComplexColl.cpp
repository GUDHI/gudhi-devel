#include <gudhi/SparseMsMatrix.h>
#include <gudhi/Fake_simplex_tree.h>

#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Random.h>
#include <cmath>

#include <iostream>
#include <string>
#include <vector>
#include <limits>  // for std::numeric_limits
#include <sstream>      // for istringstream
#include <algorithm>    // for copy, copy_n
#include <iterator>     // for istream_iterator<>, ostream_iterator<>
#include <set>

using Point = CGAL::Epick_d< CGAL::Dimension_tag<20> >::Point_d;
using Filtration_value = Fake_simplex_tree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

using Fake_simplex_tree  = Gudhi::Fake_simplex_tree ;
using Vertex             = Fake_simplex_tree::Vertex;

using Vector_of_points = std::vector<Point>;

// template<typename Set> std::set<Set> powerset(const Set& s, size_t n)
// {
//     typedef typename Set::const_iterator SetCIt;
//     typedef typename std::set<Set>::const_iterator PowerSetCIt;
//     std::set<Set> res;
//     if(n > 0) {
//         std::set<Set> ps = powerset(s, n-1);
//         for(PowerSetCIt ss = ps.begin(); ss != ps.end(); ss++)
//             for(SetCIt el = s.begin(); el != s.end(); el++) {
//                 Set subset(*ss);
//                 subset.insert(*el);
//                 res.insert(subset);
//             }
//         res.insert(ps.begin(), ps.end());
//     } else
//         res.insert(Set());
//     return res;
// }
// template<typename Set> std::set<Set> powerset(const Set& s)
// {
//     return powerset(s, s.size());
// }
 

int main()
{
	
	Vector_of_points points;
	//Vector_of_points noisePoints;
	Fake_simplex_tree stree;
	std::vector<std::vector<int>> contents;
    std::cout << "Please enter the #vertices and then the vertices of the maximal simplices in the following format" <<std::endl;
    std::cout << "3" << std::endl;
    std::cout << "1 2 3" << std::endl;
    std::cout << "Enter s to stop" << std::endl;
    int number;
    while (std::cin >> number)
    {
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip eol
        std::string line;
        std::getline(std::cin, line);
        if (std::cin)
        {
            contents.emplace_back(number);
            std::istringstream iss(line);
            std::copy_n(std::istream_iterator<double>(iss), number, contents.back().begin());
        }
        else
        {
            return 255;
        }
    }

    // if (!std::cin.eof())
    //     std::cout << "Warning: end of file not reached\n";
    
    for (auto& row : contents)
    {
        stree.insert_simplex_and_subfaces(row);
        //std::copy(row.begin(), row.end(), std::ostream_iterator<double>(std::cout," "));
        //std::cout << "\n";
    }
    // auto result = powerset(contents[0]);

    // for (auto&& subset: result) 
    // {
    //     std::cout << "{ ";
    //     char const* prefix = "";
    //     for (auto&& e: subset) {
    //         std::cout << prefix << e;
    //         prefix = ", ";
    //     }
    //     std::cout << " }\n";
    // }
	
	std::cout << "Input complex is of dimension " << stree.dimension() << " with " << stree.num_simplices() << " maximal simplices, " << stree.filtration_simplex_range().size() << " simplices and " << stree.num_vertices() << " vertices." << std::endl;

	clock_t stree_formed = clock();
	std::cout << "Now going for matrix formation" << std::endl;

	SparseMsMatrix mat(stree);
	clock_t matrix_formed = clock();
	std::cout << "Matrix formed ... Now going for collapse" << std::endl;

	Fake_simplex_tree coll_tree = mat.collapsed_tree();
	clock_t collapse_done = clock();

	std::cout << "Collapse done !" << std::endl;

	std::cout << "Time for formation of Matrix : " << (matrix_formed - stree_formed)/CLOCKS_PER_SEC << " seconds" << std::endl;
	std::cout << "Time for Collapse : " << (collapse_done - matrix_formed)/CLOCKS_PER_SEC  << " seconds" << std::endl;

	std::cout << "Collapsed Rips complex is of dimension " << coll_tree.dimension() << " with " << coll_tree.num_simplices() << " maximal simplices " << coll_tree.filtration_simplex_range().size() << " simplices and "  << coll_tree.num_vertices() << " vertices." << std::endl;


	return 0;
}


