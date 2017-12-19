#include <gudhi/SparseMsMatrix.h>
#include <gudhi/Fake_simplex_tree.h>

#include <gudhi/Alpha_complex.h>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Random.h>
#include <cmath>
#include <chrono>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>  // for std::numeric_limits

using Kernel = CGAL::Epick_d< CGAL::Dimension_tag<8> >;
using Point  = Kernel::Point_d;
using Filtration_value = Fake_simplex_tree::Filtration_value;

using Alpha_complex = Gudhi::alpha_complex::Alpha_complex<Kernel> ;

using Fake_simplex_tree = Gudhi::Fake_simplex_tree ;

using Vector_of_points = std::vector<Point>;

void generate_points_sphere(Vector_of_points& W, int nbP, int dim)
{
	CGAL::Random_points_on_sphere_d<Point> rp(dim, 4);
	for (int i = 0; i < nbP; i++)
		W.push_back(*rp++);
}
void generate_points_ball(Vector_of_points& W, int nbP, int dim, int radius)
{
	CGAL::Random_points_in_ball_d<Point> rp(dim, radius); 
	for (int i = 0; i < nbP; i++)
		W.push_back(*rp++);
}

void generate_points_cube(Vector_of_points& W, int nbP, int dim)
{
	CGAL::Random_points_in_cube_d<Point> rp(dim, 6);
	for (int i = 0; i < nbP; i++)
		W.push_back(*rp++);
}



int main()
{
	
	Vector_of_points points;
	
	int originalDim, collDim;
	long originalNumSimp, colNumSimp;
	int originalNumVert, colNumVert;
	long originalNumMxSimp, colNumMxSimp;
	
	int num_pts = 300;
	int dimBall	= 2;
	int radius  = 3;
	double threshold =  0.03;

	//generate_points_sphere(points,num_pts,3);
	generate_points_ball(points,num_pts,dimBall,radius); 
	//add_point_vectors(points, noisePoints, num_pts);
	
	std::cout << "Number of points : " << num_pts << " and " << "Threshold value : " << threshold << std::endl;
	Alpha_complex alpha_complex_from_points(points);
	
	// std::ofstream myfileOrig ("./output/rips1_original.xls");
	// std::ofstream myfileCol  ("./output/rips1_collapsed.xls");
	std::ofstream myfileDetl ("./output/AlphaComplexOn2-Ball(3)AmbDim8.xls");

	// myfileOrig << "Threshold, Dimension, NumOfSimplices" << std::endl;
	// myfileOrig << "Threshold, Dimension, NumOfSimplices"<< std::endl;
	myfileDetl << num_pts << " points chosen randomly from "<< dimBall <<"-ball of radius " << radius << std::endl;
	myfileDetl << "Thresold, InitialNumVertices, AfterCollNumVert, InitialDimension, CollDimension, InitialNumSimp, CollapsedNumSimplices, TimeInMS" << std::endl;

	while(threshold <= 1.1)
	{
		Fake_simplex_tree stree;
		if(alpha_complex_from_points.create_complex(stree, threshold))
		{
			// auto stree_formed  = std::chrono::high_resolution_clock::now();
			std::cout << "Simplex tree created... Next stop matrix formation" << std::endl;

			SparseMsMatrix mat(stree);
			auto matrix_formed  = std::chrono::high_resolution_clock::now();
			std::cout << "Matrix formed... Next action COLLAPSE!!" << std::endl;

			Fake_simplex_tree coll_tree = mat.collapsed_tree();
			auto collapse_done = std::chrono::high_resolution_clock::now();

			std::cout << "Collapse done !" << std::endl;
			auto collapseTime = std::chrono::duration<double, std::milli>(collapse_done- matrix_formed).count();
			// std::cout << "Time for formation of Matrix : " << (matrix_formed - stree_formed)/CLOCKS_PER_SEC << " seconds" << std::endl;
			std::cout << "Time for Collapse : " << collapseTime << " ms\n" << std::endl;
			
			originalDim = stree.dimension();
			collDim 	= coll_tree.dimension();
			
			originalNumVert = stree.num_vertices();
			colNumVert		= coll_tree.num_vertices();

			originalNumMxSimp 	= stree.num_simplices();
			colNumMxSimp 		= coll_tree.num_simplices();

			originalNumSimp = stree.filtration_simplex_range().size();
			colNumSimp		= coll_tree.filtration_simplex_range().size();
			

			std::cout << "Rips complex is of dimension " << originalDim << " with " << originalNumMxSimp << " maximal simplices, " << originalNumSimp << " simplices and " << originalNumVert << " vertices." << std::endl;
			std::cout << "Collapsed Rips complex is of dimension " << collDim << " with " <<  colNumMxSimp << " maximal simplices, " <<  colNumSimp << " simplices and " << colNumVert << " vertices." << std::endl;

			std::cout << "** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** " << std::endl;


			// myfileOrig << threshold << "," << originalDim << "," << originalNumSimp << std::endl;
			// myfileOrig << threshold << "," << collDim << "," << colNumMxSimp << std::endl;
			myfileDetl << threshold << "," << originalNumVert << "," << colNumVert << "," << originalDim << "," << collDim << "," << originalNumSimp  <<  "," << colNumSimp << "," << collapseTime << std::endl;

			threshold = threshold+0.01;
		}	
		
	}
	// myfileOrig.close();
	// myfileCol.close();
	myfileDetl.close();

	return 0;
}
