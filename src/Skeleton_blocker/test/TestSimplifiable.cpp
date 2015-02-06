 /*    This file is part of the Gudhi Library. The Gudhi library 
  *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
  *    library for computational topology.
  *
  *    Author(s):       David Salinas
  *
  *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
  *
  *    This program is free software: you can redistribute it and/or modify
  *    it under the terms of the GNU General Public License as published by
  *    the Free Software Foundation, either version 3 of the License, or
  *    (at your option) any later version.
  *
  *    This program is distributed in the hope that it will be useful,
  *    but WITHOUT ANY WARRANTY; without even the implied warranty of
  *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *    GNU General Public License for more details.
  *
  *    You should have received a copy of the GNU General Public License
  *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
  */


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include "gudhi/Test.h"
//#include "Skeleton_blocker/Simplex.h"
#include "gudhi/Skeleton_blocker.h"


using namespace std;

using namespace Gudhi;

using namespace skbl;

template<typename ComplexType> class Skeleton_blocker_sub_complex;
typedef Skeleton_blocker_complex<Skeleton_blocker_simple_traits> Complex;
typedef Complex::Vertex_handle Vertex_handle;
typedef Complex::Root_vertex_handle Root_vertex_handle;
typedef Skeleton_blocker_simplex<Vertex_handle> Simplex_handle;
// true iff v \in complex
bool assert_vertex(Complex &complex,Vertex_handle v){
	Simplex_handle simplex(v);
	bool test = complex.contains(simplex);
	assert(test);
	return test;
}

// true iff the blocker (a,b,c) is in complex
bool assert_blocker(Complex &complex,Root_vertex_handle a,Root_vertex_handle b,Root_vertex_handle c){

	return complex.contains_blocker(Simplex_handle(*complex.get_address(a),*complex.get_address(b),*complex.get_address(c)));
	//return complex.contains_blocker((a),(b),(c));
}

void build_complete(int n,Complex& complex){
	complex.clear();
	for(int i=0;i<n;i++)
		complex.add_vertex();
	for(int i=0;i<n;i++)
		for(int j=0;j<i;j++)
			complex.add_edge(Vertex_handle(i),Vertex_handle(j));
}

bool test_contraction1(){
	enum { a, b, x, y, z, n };
	Complex complex(n);
	build_complete(n,complex);
	complex.remove_edge(static_cast<Vertex_handle>(b), static_cast<Vertex_handle>(z));
	complex.add_blocker(Simplex_handle(static_cast<Vertex_handle>(a), static_cast<Vertex_handle>(x),
	                                   static_cast<Vertex_handle>(y)));
	complex.add_blocker(Simplex_handle(static_cast<Vertex_handle>(b), static_cast<Vertex_handle>(x),
	                                   static_cast<Vertex_handle>(y)));

	// Print result
	cerr << "complex before complex"<< complex.to_string()<<endl;

	cerr <<endl<<endl;
	complex.contract_edge(static_cast<Vertex_handle>(a),static_cast<Vertex_handle>(b));
	// Print result
	cerr << "ContractEdge(0,1)\n";
	PRINT(complex.to_string());

	// verification
	for (int i=0;i<5;i++)
		if (i!=1) assert_vertex(complex, static_cast<Vertex_handle>(i));
	bool test1 = !complex.contains_edge(static_cast<Vertex_handle>(a),static_cast<Vertex_handle>(b));
	bool test2 = assert_blocker(complex,Root_vertex_handle(a),Root_vertex_handle(x),Root_vertex_handle(y));
	bool test3 = complex.num_edges()==6;
	bool test4 = complex.num_blockers()==1;
	Simplex_handle sigma;
	sigma.add_vertex(static_cast<Vertex_handle>(a));
	sigma.add_vertex(static_cast<Vertex_handle>(x));
	sigma.add_vertex(static_cast<Vertex_handle>(y));
	sigma.add_vertex(static_cast<Vertex_handle>(z));
	bool test5 = !(complex.contains(sigma));
	return test1&&test2&&test3&&test4&&test5;
}


bool test_contraction2(){
	enum { a, b, x, y, z, n };
	Complex complex(n);
	build_complete(n,complex);
	complex.remove_edge(static_cast<Vertex_handle>(b),static_cast<Vertex_handle>(x));
	Simplex_handle blocker;
	blocker.add_vertex(static_cast<Vertex_handle>(a));
	blocker.add_vertex(static_cast<Vertex_handle>(y));
	blocker.add_vertex(static_cast<Vertex_handle>(z));

	complex.add_blocker(blocker);

	// Print result
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;
	complex.contract_edge(static_cast<Vertex_handle>(a),static_cast<Vertex_handle>(b));

	cerr << "complex.ContractEdge(a,b)"<< complex.to_string();

	cerr <<endl<<endl;

	// there should be one blocker (a,c,d,e) in the complex
	bool test ;
	test = complex.contains_blocker(Simplex_handle(static_cast<Vertex_handle>(a), static_cast<Vertex_handle>(x),
	                                               static_cast<Vertex_handle>(y),static_cast<Vertex_handle>(z)));
	test = test && complex.num_blockers()==1;
	return test;
}

bool test_link_condition1(){

	Complex complex(0);
	// Build the complexes
	build_complete(4,complex);
	complex.add_blocker(Simplex_handle(static_cast<Vertex_handle>(0), static_cast<Vertex_handle>(1), static_cast<Vertex_handle>(2)));


	// Print result
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;

	bool weak_link_condition = complex.link_condition(Vertex_handle(1),Vertex_handle(2),true);

	bool strong_link_condition = complex.link_condition(Vertex_handle(1),Vertex_handle(2),false);

	return weak_link_condition && !strong_link_condition;
}


bool test_collapse0(){
	Complex complex(5);
	build_complete(4,complex);
	complex.add_vertex();
	complex.add_edge(static_cast<Vertex_handle>(2), static_cast<Vertex_handle>(4));
	complex.add_edge(static_cast<Vertex_handle>(3), static_cast<Vertex_handle>(4));
	// Print result
	cerr << "initial complex :\n"<< complex.to_string();
	cerr <<endl<<endl;

	Simplex_handle simplex_123(static_cast<Vertex_handle>(1), static_cast<Vertex_handle>(2), static_cast<Vertex_handle>(3));
	complex.remove_star(simplex_123);
	cerr << "complex.remove_star(1,2,3):\n"<< complex.to_string();
	cerr <<endl<<endl;

	// verification
	bool blocker123_here = complex.contains_blocker(simplex_123);
	cerr <<"----> Ocomplex \n";
	return blocker123_here;
}


bool test_collapse1(){
	Complex complex(5);
	build_complete(4,complex);
	complex.add_blocker(Simplex_handle(Vertex_handle(0),Vertex_handle(1),Vertex_handle(2),Vertex_handle(3)));
	// Print result
	cerr << "initial complex :\n"<< complex.to_string();
	cerr <<endl<<endl;

	Simplex_handle simplex_123(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3));
	complex.remove_star(simplex_123);
	cerr << "complex.remove_star(1,2,3):\n"<< complex.to_string();
	cerr <<endl<<endl;

	// verification
	bool res = complex.contains_blocker(simplex_123);
	res = res && complex.num_blockers()==1;
	cerr <<"----> Ocomplex \n";
	return res;
}

bool test_collapse2(){
	Complex complex(5);
	build_complete(4,complex);
	complex.add_vertex();
	complex.add_edge(Vertex_handle(1),Vertex_handle(4));
	complex.add_edge(Vertex_handle(2),Vertex_handle(4));
	complex.add_edge(Vertex_handle(3),Vertex_handle(4));
	complex.add_blocker(Simplex_handle(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3),Vertex_handle(4)));
	// Print result
	cerr << "initial complex :\n"<< complex.to_string();
	cerr <<endl<<endl;

	Simplex_handle sigma(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3));
	complex.remove_star(sigma);
	cerr << "complex.remove_star(1,2,3):\n"<< complex.to_string();
	cerr <<endl<<endl;

	// verification
	bool blocker_removed = !complex.contains_blocker(Simplex_handle(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3),Vertex_handle(4)));
	bool blocker123_here = complex.contains_blocker(sigma);
	return blocker_removed && blocker123_here;
}

bool test_collapse3(){
	Complex complex(5);
	build_complete(4,complex);
	complex.add_vertex();
	complex.add_edge(Vertex_handle(1),Vertex_handle(4));
	complex.add_edge(Vertex_handle(2),Vertex_handle(4));
	complex.add_edge(Vertex_handle(3),Vertex_handle(4));
	complex.add_blocker(Simplex_handle(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3),Vertex_handle(4)));
	// Print result
	cerr << "initial complex:\n"<< complex.to_string();
	cerr <<endl<<endl;

	complex.remove_star(static_cast<Vertex_handle>(2));
	cerr << "complex after remove star of 2:\n"<< complex.to_string();

	bool blocker134_here = complex.contains_blocker(Simplex_handle(Vertex_handle(1),Vertex_handle(3),Vertex_handle(4)));
	bool blocker1234_here = complex.contains_blocker(Simplex_handle(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3),Vertex_handle(4)));
	return blocker134_here && !blocker1234_here;
}

bool test_add_simplex(){
	Complex complex(5);
	build_complete(4,complex);
	complex.add_blocker(Simplex_handle(Vertex_handle(0),Vertex_handle(1),Vertex_handle(2)));
	// Print result
	cerr << "initial complex:\n"<< complex.to_string();
	cerr <<endl<<endl;

	complex.add_simplex(Simplex_handle(Vertex_handle(0),Vertex_handle(1),Vertex_handle(2),Vertex_handle(3)));
	cerr << "complex after add_simplex:\n"<< complex.to_string();

	return complex.num_blockers()==0;
}

bool test_add_simplex2(){
	Complex complex;
	build_complete(4,complex);
	// Print result
	cerr << "initial complex:\n"<< complex.to_string();
	cerr <<endl<<endl;

	Complex copy(complex.num_vertices());

	std::vector<Simplex_handle> simplices(complex.simplex_range().begin(),complex.simplex_range().end());
	sort(simplices.begin(),simplices.end(),[&](const Simplex_handle& s1,const Simplex_handle& s2){
		return s1.dimension()<s2.dimension();
	});
	for(const auto & simplex : simplices){
		if(!copy.contains(simplex) && simplex.dimension()==1)
			copy.add_edge(simplex.first_vertex(),simplex.last_vertex());
		if(!copy.contains(simplex) && simplex.dimension()>1)
			copy.add_simplex(simplex);
	}


	cerr << "complex after add_simplex:\n"<< copy.to_string();


	return complex.num_blockers()==copy.num_blockers() &&
			complex.num_edges()==copy.num_edges() &&
			complex.num_vertices()==copy.num_vertices();
}



bool test_remove_popable_blockers(){
	Complex complex;
	build_complete(4,complex);
	complex.add_vertex();
	complex.add_edge(Vertex_handle(3),Vertex_handle(4));
	complex.add_edge(Vertex_handle(2),Vertex_handle(4));
	Simplex_handle sigma1=Simplex_handle(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3));
	Simplex_handle sigma2=Simplex_handle(Vertex_handle(2),Vertex_handle(3),Vertex_handle(4));

	complex.add_blocker(sigma1);
	complex.add_blocker(sigma2);
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;
	cerr << "complex.RemovePopableBlockers();"<<endl;
	complex.remove_popable_blockers();
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;

	bool test1 = (complex.num_blockers()==1);


	// test 2
	complex.clear();
	build_complete(4,complex);
	complex.add_vertex();
	complex.add_vertex();
	complex.add_edge(Vertex_handle(3),Vertex_handle(4));
	complex.add_edge(Vertex_handle(2),Vertex_handle(4));
	complex.add_edge(Vertex_handle(2),Vertex_handle(5));
	complex.add_edge(Vertex_handle(3),Vertex_handle(5));
	complex.add_edge(Vertex_handle(4),Vertex_handle(5));
	sigma1=Simplex_handle(Vertex_handle(1),Vertex_handle(2),Vertex_handle(3));
	sigma2=Simplex_handle(Vertex_handle(2),Vertex_handle(3),Vertex_handle(4));

	complex.add_blocker(sigma1);
	complex.add_blocker(sigma2);
	cerr << "complex complex"<< complex.to_string();
	cerr <<endl<<endl;	cerr << "complex.RemovePopableBlockers();"<<endl;
	complex.remove_popable_blockers();
	cerr << "complex complex"<< complex.to_string();

	cerr <<endl<<endl;
	bool test2 = (complex.num_blockers()==0);
	return test1&&test2;
}



int main (int argc, char *argv[])
{
	Tests tests_simplifiable_complex;
	tests_simplifiable_complex.add("Test contraction 1",test_contraction1);
	tests_simplifiable_complex.add("Test contraction 2",test_contraction2);
	tests_simplifiable_complex.add("Test Link condition 1",test_link_condition1);
	tests_simplifiable_complex.add("Test remove popable blockers",test_remove_popable_blockers);


	tests_simplifiable_complex.add("Test collapse 0",test_collapse0);
	tests_simplifiable_complex.add("Test collapse 1",test_collapse1);
	tests_simplifiable_complex.add("Test collapse 2",test_collapse2);
	tests_simplifiable_complex.add("Test collapse 3",test_collapse3);

	tests_simplifiable_complex.add("Test add simplex",test_add_simplex);
	tests_simplifiable_complex.add("Test add simplex2",test_add_simplex2);


	tests_simplifiable_complex.run();

	if(tests_simplifiable_complex.run())
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}
