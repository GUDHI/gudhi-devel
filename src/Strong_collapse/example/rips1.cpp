#include <gudhi/MsMatrix.h>


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

using Point = CGAL::Epick_d< CGAL::Dimension_tag<10> >::Point_d;
using Filtration_value = Simplex_tree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

using Vector_of_points = std::vector<Point>;

void generate_points_sphere(Vector_of_points& W, int nbP, int dim)
{
	CGAL::Random_points_on_sphere_d<Point> rp(dim, 1);
	for (int i = 0; i < nbP; i++)
		W.push_back(*rp++);
}
void generate_points_ball(Vector_of_points& W, int nbP, int dim)
{
	CGAL::Random_points_in_ball_d<Point> rp(dim, 1);
	for (int i = 0; i < nbP; i++)
		W.push_back(*rp++);
}



int main()
{
	Vector_of_points points;
	int num_pts = 15;
	generate_points_ball(points,num_pts,4);

	double threshold = 0.8;
	std::cout << "Number of points : " << num_pts << " and " << "Threshold value : " << threshold << std::endl;
	Rips_complex rips_complex_from_points(points, threshold, Gudhi::Euclidean_distance());

	for(int dime = 3 ; dime <= 6 ; ++dime)
	{
		Simplex_tree stree;
		rips_complex_from_points.create_complex(stree, dime);

		std::cout << "Rips complex is of dimension " << dime << " with " << stree.num_simplices() << " simplices and " << stree.num_vertices() << " vertices." << std::endl;

		clock_t stree_formed = clock();
		std::cout << "Simplex tree formed ... Now going for matrix formation" << std::endl;

		MsMatrix mat(stree);
		clock_t matrix_formed = clock();
		std::cout << "Matrix formed ... Now going for collapse" << std::endl;

		Simplex_tree coll_tree = mat.collapsed_tree();
		clock_t collapse_done = clock();

		std::cout << "Collapse done !" << std::endl;

		std::cout << "Time for formation of Matrix : " << (matrix_formed - stree_formed)/CLOCKS_PER_SEC << " seconds" << std::endl;
		std::cout << "Time for Collapse : " << (collapse_done - matrix_formed)/CLOCKS_PER_SEC  << " seconds" << std::endl;

		std::cout << "Collapsed Rips complex is of " << coll_tree.num_simplices() << " simplices and " << coll_tree.num_vertices() << " vertices." << std::endl;

		std::cout << "** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** " << std::endl;

		FILE * myfile;

		myfile = fopen("rips1_original.txt","a");
		fprintf (myfile, "%d %lu\n" , dime , stree.num_simplices() );
		fclose(myfile);

		myfile = fopen("rips1_collapsed.txt","a");
		fprintf (myfile, "%d %lu\n" , dime , coll_tree.num_simplices() );
		fclose(myfile);

		myfile = fopen("rips1_details.txt","a");
		fprintf (myfile, "%d %lu %lu %lu\n" , dime , stree.num_simplices() , coll_tree.num_simplices() , (collapse_done - matrix_formed)/CLOCKS_PER_SEC );
		fclose(myfile);
	}

	return 0;
}
