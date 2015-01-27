/*
 * Model.h
 *
 *  Created on: Aug 25, 2014
 *      Author: david
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <fstream>
#include <limits>
#include "gudhi/Clock.h"
#include "utils/UI_utils.h"
#include "utils/Lloyd_builder.h"
#include "utils/Rips_builder.h"
#include "utils/K_nearest_builder.h"
#include "utils/Vertex_collapsor.h"
#include "utils/Edge_collapsor.h"
#include "utils/Edge_contractor.h"
#include "utils/Persistence_compute.h"
#include "utils/Critical_points.h"

#include "gudhi/Skeleton_blocker/Skeleton_blocker_simple_geometric_traits.h"
#include "gudhi/Skeleton_blocker_geometric_complex.h"

#include "gudhi/Off_reader.h"

#include "Complex_typedefs.h"


#include <CGAL/Euclidean_distance.h>


template<typename Complex>
class CGAL_geometric_flag_complex_wrapper{
	Complex& complex_;
	typedef typename Complex::Vertex_handle Vertex_handle;
	typedef typename Complex::Point Point;

	const bool load_only_points_;

public:
	CGAL_geometric_flag_complex_wrapper(Complex& complex,bool load_only_points = false):
		complex_(complex),
		load_only_points_(load_only_points){
	}

	void init(int dim,int num_vertices,int num_max_faces,int num_edges) const{
	}

	void point(const std::vector<double>& coords){
		Point p(coords.size(),coords.begin(),coords.end());
		complex_.add_vertex(p);
	}

	void maximal_face(std::vector<int> vertices){
		if (!load_only_points_){
			for (int i = 0; i<vertices.size()-1 ; ++i)
				for (int j = i+1; j<vertices.size()-1 ; ++j)
					complex_.add_edge(Vertex_handle(vertices[i]),Vertex_handle(vertices[j]));
		}
	}
	void done() const{}
};


class Model{

public:
	Complex complex_;
	typedef Complex::Vertex_handle Vertex_handle;
	typedef Complex::CVI CVI;


	Model():complex_(){
	}

public:
	void off_file_open(const std::string& name_file){
		UIDBGMSG("load off file",name_file);
		complex_.clear();
		CGAL_geometric_flag_complex_wrapper<Complex> read_wraper(complex_);
		Gudhi::read_off(name_file,read_wraper);
	}

	void off_points_open(const std::string& name_file){
		UIDBGMSG("load off points",name_file);
		complex_.clear();
		CGAL_geometric_flag_complex_wrapper<Complex> read_wraper(complex_);
		Gudhi::read_off(name_file,read_wraper);
	}

	void off_file_save(const std::string& name_file){
		UIDBG("save off file");
	}

	void off_points_save(const std::string& name_file){
		UIDBG("save off off_points_save");
	}

	// point sets operations
	void uniform_noise(double amplitude){
		UIDBG("unif noise");
		for (auto v : complex_.vertex_range())
			complex_.point(v) = add_uniform_noise(complex_.point(v),amplitude);
	}

private:
	Point add_uniform_noise(const Point& point,double amplitude){
		std::vector<double> new_point(point.dimension());
		for(int i  = 0 ; i < point.dimension();++i){
			new_point[i] = point[i] + (rand() % 2 - .5) * amplitude;
		}
		return Point(point.dimension(), new_point.begin(),new_point.end());
	}

public:

	void lloyd(int num_iterations,int num_closest_neighbors){
		UIDBG("lloyd");
		Lloyd_builder<Complex> lloyd_builder(complex_,1);
	}

	double squared_eucl_distance(const Point& p1,const Point& p2) const{
		return Geometry_trait::Squared_distance_d()(p1,p2);
	}

	// complex operations from points
	void build_rips(double alpha){
		UIDBG("build_rips");
		Rips_builder<Complex> rips_builder(complex_,alpha);
	}

	void build_k_nearest_neighbors(unsigned k){
		UIDBG("build_k_nearest");
		complex_.keep_only_vertices();
		K_nearest_builder<Complex> k_nearest_builder(complex_,k);
	}

	void build_delaunay(){
		UIDBG("build_delaunay");
		complex_.keep_only_vertices();
	}



	void contract_edges(unsigned num_contractions){
		Clock c;
		Edge_contractor<Complex> contractor(complex_,num_contractions);
		std::cout<<"Time to simplify: "<<c.num_seconds()<<"s"<<std::endl;
	}


	void collapse_vertices(unsigned num_collapses){

		auto old_num_vertices = complex_.num_vertices();
		Vertex_collapsor<Complex> collapsor(complex_,complex_.num_vertices());
		UIDBGMSG("num vertices collapsed:",old_num_vertices - complex_.num_vertices());
	}

	void collapse_edges(unsigned num_collapses){
		Edge_collapsor<Complex> collapsor(complex_,complex_.num_edges());
	}


	void show_graph_stats(){
		std::cout << "++++++ Graph stats +++++++"<< std::endl;
		std::cout << "Num vertices : " << complex_.num_vertices()<<std::endl;
		std::cout << "Num edges : " << complex_.num_edges()<<std::endl;
		std::cout << "Num connected components : " << complex_.num_connected_components()<<std::endl;
		std::cout << "Min/avg/max degree : " << min_degree()<<"/"<<avg_degree()<<"/"<<max_degree()<<std::endl;
		std::cout << "Num connected components : " << complex_.num_connected_components()<<std::endl;
		std::cout << "Num connected components : " << complex_.num_connected_components()<<std::endl;
		std::cout << "+++++++++++++++++++++++++"<< std::endl;
	}

private:
	int min_degree() const{
		int res = (std::numeric_limits<int>::max)();
		for(auto v : complex_.vertex_range())
			res= (std::min)(res,complex_.degree(v));
		return res;
	}

	int max_degree() const{
		int res = 0;
		for(auto v : complex_.vertex_range())
			res= (std::max)(res,complex_.degree(v));
		return res;
	}

	int avg_degree() const{
		int res = 0;
		for(auto v : complex_.vertex_range())
			res+= complex_.degree(v);
		return res / complex_.num_vertices();
	}

public:



	void show_complex_stats(){
		std::cout << "++++++ Mesh stats +++++++"<< std::endl;
		std::cout << "Num vertices : " << complex_.num_vertices()<<std::endl;
		std::cout << "Num edges : " << complex_.num_edges()<<std::endl;
		std::cout << "Num connected components : " << complex_.num_connected_components()<<std::endl;
		std::cout << "+++++++++++++++++++++++++"<< std::endl;

	}

	void show_complex_dimension(){
		unsigned num_simplices = 0;
		int euler = 0;
		int dimension = 0;
		Clock clock;
		for(const auto &s : complex_.simplex_range()){
			num_simplices++;
			dimension = (std::max)(s.dimension(),dimension);
			if(s.dimension()%2==0)
				euler+=1;
			else
				euler-=1;
		}
		clock.end();
		std::cout << "++++++ Mesh dimension +++++++"<< std::endl;
		std::cout << "Dimension : " << dimension<<std::endl;
		std::cout << "Euler characteristic : " << euler<<std::endl;
		std::cout << "Num simplices : " << num_simplices <<std::endl;
		std::cout << "Total time: " << clock <<std::endl;
		std::cout << "Time per simplex: " << clock.num_seconds()/num_simplices <<" s"<<std::endl;
		std::cout << "+++++++++++++++++++++++++"<< std::endl;
	}


	void show_homology_group(){
#ifdef _WIN32
		std::cout << "Works only on linux x64 for the moment\n";
#else
		Clock clock;
		run_chomp();
		clock.end();
#endif
	}

	void show_euler_characteristic(){
		unsigned num_simplices = 0;
		int euler = 0;
		int dimension = 0;
		for(const auto &s : complex_.simplex_range()){
			num_simplices++;
			dimension = (std::max)(s.dimension(),dimension);
			if(s.dimension()%2==0)
				euler+=1;
			else
				euler-=1;
		}
		std::cout << "Saw "<<num_simplices<<" simplices with maximum dimension " << dimension<<std::endl;
		std::cout << "The euler characteristic is : " << euler<<std::endl;
	}

	void show_persistence(	int p,double threshold,int max_dim,double min_pers){
		Persistence_compute<Complex> persistence(complex_,std::cout,Persistence_params(p,threshold,max_dim,min_pers));
	}

	void show_critical_points(double max_distance){
		Critical_points<Complex> critical_points(complex_,std::cout,max_distance);
	}


private:
	void run_chomp(){
		save_complex_in_file_for_chomp();
		std::cout << "Call CHOMP library\n";
		system("../src/utils/homsimpl chomp.sim");
	}

	void save_complex_in_file_for_chomp(){
		std::ofstream file;
		file.open("chomp.sim");
		for(const auto &s : complex_.simplex_range()){
			bool first = true;
			file<<"(";
			for(auto x : s){
				if(first) first = false;
				else file<<",";
				file << x;
			}
			file<<")\n";
		}
	}
public:


	unsigned num_vertices() const{
		return  complex_.num_vertices();
	}
};

#endif /* MODEL_H_ */
