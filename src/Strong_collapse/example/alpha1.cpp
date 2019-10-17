#include <gudhi/MsMatrix.h>

#include <gudhi/Alpha_complex.h>
#include <CGAL/Epick_d.h>

#include <CGAL/point_generators_d.h>
#include <CGAL/Random.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>

using Kernel = CGAL::Epick_d< CGAL::Dimension_tag<3> >;
using Point = Kernel::Point_d;
using Vector_of_points = std::vector<Point>;

void generate_points_sphere(Vector_of_points& W, int nbP, int dim)
{
	CGAL::Random_points_on_sphere_d<Point> rp(dim, 1);
	for (int i = 0; i < nbP; i++)
		W.push_back(*rp++);
}

void usage(int nbArgs, char * const progName)
{
	std::cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
	std::cerr << "Usage: " << progName << " [alpha_square_max_value]\n";
	std::cerr << "       i.e.: " << progName << " 60.0\n";
	exit(-1);  // ----- >>
}

int main(int argc, char **argv)
{
	if ((argc != 1) && (argc != 2)) usage(argc, (argv[0] - 1));

	// Delaunay complex if alpha_square_max_value is not given by the user.
	double alpha_square_max_value {std::numeric_limits<double>::infinity()};

	if (argc == 2)
		alpha_square_max_value = atof(argv[1]);

	int num_pts;
	Vector_of_points points;

	for(num_pts = 50 ; num_pts <= 5000 ; num_pts += 50)
	{
		points.clear();
		generate_points_sphere(points,num_pts,2);

		Gudhi::alpha_complex::Alpha_complex<Kernel> alpha_complex_from_points(points);

		Gudhi::Simplex_tree<> simplex;
		if (alpha_complex_from_points.create_complex(simplex, alpha_square_max_value)) 
		{
			std::cout << "Original Alpha complex is of " << simplex.num_simplices() << " simplices - " << simplex.num_vertices() << " vertices." << std::endl;

			clock_t stree_formed = clock();
		 	std::cout << "simplex tree formed ... Now going for matrix formation" << std::endl;

		 	MsMatrix mat(simplex);
		 	clock_t matrix_formed = clock();
		 	std::cout << "Matrix formed ... Now going for collapse" << std::endl;

		 	Simplex_tree coll_tree = mat.collapsed_tree();
		 	clock_t collapse_done = clock();

		 	std::cout << "Collapse done !" << std::endl;

			std::cout << "Time for formation of Matrix : " << (matrix_formed - stree_formed)/CLOCKS_PER_SEC << " seconds" << std::endl;
			std::cout << "Time for Collapse : " << (collapse_done - matrix_formed)/CLOCKS_PER_SEC  << " seconds" << std::endl;

			std::cout << "Collapsed Alpha complex is of " << coll_tree.num_simplices() << " simplices and " << coll_tree.num_vertices() << " vertices." << std::endl;

			std::cout << "** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** " << std::endl;

			FILE * myfile;

			myfile = fopen("alpha1_original.txt","a");
			fprintf (myfile, "%d %lu\n" , num_pts , simplex.num_simplices() );
			fclose(myfile);

			myfile = fopen("alpha1_collapsed.txt","a");
			fprintf (myfile, "%d %lu\n" , num_pts , coll_tree.num_simplices() );
			fclose(myfile);

			myfile = fopen("alpha1_details.txt","a");
			fprintf (myfile, "%d %lu %lu %lu\n" , num_pts , simplex.num_simplices() , coll_tree.num_simplices() , (collapse_done - matrix_formed)/CLOCKS_PER_SEC);
			fclose(myfile);
		}
		else
		{
			std::cout << "Create_complex failed !" << std::endl;
		}
	}

	return 0;
}
