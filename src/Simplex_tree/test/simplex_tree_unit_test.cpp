#define BOOST_TEST_MODULE const_string test
#include <boost/test/included/unit_test.hpp>
#include <boost/system/error_code.hpp>
#include <boost/chrono/thread_clock.hpp>
#include <iostream>
#include <string>

#include <utility> // std::pair, std::make_pair

#include <cmath> // float comparison
#include <limits>

#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/reader_utils.h"
#include "gudhi/Simplex_tree.h"


using namespace Gudhi;

typedef Simplex_tree<> typeST;
typedef std::pair< typeST::Simplex_handle, bool > typePairSimplexBool;
typedef std::vector< Vertex_handle > typeVectorVertex;
typedef std::pair<typeVectorVertex, Filtration_value> typeSimplex;

const Vertex_handle DEFAULT_VERTEX_HANDLE = (const Vertex_handle)-1;

void test_empty_simplex_tree(typeST& tst)
{
	BOOST_CHECK( tst.null_vertex() == DEFAULT_VERTEX_HANDLE );
	BOOST_CHECK( tst.filtration() == Filtration_value(0) );
	BOOST_CHECK( tst.num_vertices() == (size_t)0 );
	BOOST_CHECK( tst.num_simplices() == (size_t)0 );
	typeST::Siblings* STRoot = tst.root();
	BOOST_CHECK( STRoot != NULL );
	BOOST_CHECK( STRoot->oncles() == NULL );
	BOOST_CHECK( STRoot->parent() == DEFAULT_VERTEX_HANDLE );
	BOOST_CHECK( tst.dimension() == -1 );
}

void test_iterators_on_empty_simplex_tree(typeST& tst)
{
	std::cout << "Iterator on vertices: " << std::endl;
	for( auto vertex : tst.complex_vertex_range() )
	{
		std::cout << "vertice:" << vertex << std::endl;
		BOOST_CHECK( false); // shall be empty
	}
	std::cout << "Iterator on simplices: " << std::endl;
	for( auto simplex : tst.complex_simplex_range() )
	{
		BOOST_CHECK( false); // shall be empty
	}

	std::cout << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
	for( auto f_simplex : tst.filtration_simplex_range() )
	{
		BOOST_CHECK( false); // shall be empty
		std::cout << "test_iterators_on_empty_simplex_tree - filtration=" << tst.filtration(f_simplex) << std::endl;
	}
}

BOOST_AUTO_TEST_CASE( simplex_tree_when_empty )
{
	const Filtration_value DEFAULT_FILTRATION_VALUE = 0;

	// TEST OF DEFAULT CONSTRUCTOR
	std::cout << "********************************************************************" << std::endl;
	std::cout << "TEST OF DEFAULT CONSTRUCTOR" << std::endl;
	typeST st;

	test_empty_simplex_tree(st);

	test_iterators_on_empty_simplex_tree(st);
	// TEST OF EMPTY INSERTION
	std::cout << "TEST OF EMPTY INSERTION" << std::endl;
	typeVectorVertex simplexVectorEmpty;
	BOOST_CHECK(simplexVectorEmpty.empty() == true);
	typePairSimplexBool returnEmptyValue = st.insert ( simplexVectorEmpty, DEFAULT_FILTRATION_VALUE );
	BOOST_CHECK(returnEmptyValue.first == typeST::Simplex_handle(NULL));
	BOOST_CHECK(returnEmptyValue.second == true);

	test_empty_simplex_tree(st);

	test_iterators_on_empty_simplex_tree(st);
}

bool AreAlmostTheSame(float a, float b)
{
	return std::fabs(a - b) < std::numeric_limits<float>::epsilon();
}

BOOST_AUTO_TEST_CASE( simplex_tree_from_file )
{
	// TEST OF INSERTION
	std::cout << "********************************************************************" << std::endl;
	std::cout << "TEST OF SIMPLEX TREE FROM A FILE" << std::endl;
	typeST st;

	std::string inputFile("simplex_tree_for_unit_test.txt");
	std::ifstream simplex_tree_stream(inputFile.c_str());
	simplex_tree_stream >> st;

	// Display the Simplex_tree
	std::cout << "The complex contains " << st.num_simplices() << " simplices" << std::endl;
	std::cout << "   - dimension " << st.dimension() << "   - filtration " << st.filtration() << std::endl;

	// Check
	BOOST_CHECK(st.num_simplices() == 143353);
	BOOST_CHECK(st.dimension() == 3);
	BOOST_CHECK(st.filtration() == 0.4);

	int previous_size=0;
	for( auto f_simplex : st.filtration_simplex_range() )
	{
		// Size of simplex
		int size=0;
		for( auto vertex : st.simplex_vertex_range(f_simplex) )
		{
			size++;
		}
		BOOST_CHECK(AreAlmostTheSame(st.filtration(f_simplex),(0.1* size))); // Specific test: filtration = 0.1 * simplex_size
		BOOST_CHECK(previous_size <= size); // Check list is sorted (because of sorted filtrations in simplex_tree.txt)
		previous_size = size;
	}
	simplex_tree_stream.close();
}

void test_simplex_tree_contains(typeST& simplexTree, typeSimplex& simplex, int pos)
{
	auto f_simplex = simplexTree.filtration_simplex_range().begin() + pos;

	std::cout << "test_simplex_tree_contains - filtration=" << simplexTree.filtration(*f_simplex) << "||" << simplex.second << std::endl;
	BOOST_CHECK( AreAlmostTheSame(simplexTree.filtration(*f_simplex),simplex.second) );

	//typeVectorVertex::iterator simplexIter = simplex.first.end()-1;
	int simplexIndex=simplex.first.size()-1;
	for( auto vertex : simplexTree.simplex_vertex_range(*f_simplex) )
	{
		//std::cout << "test_simplex_tree_contains - vertex=" << vertex << "||" << *simplexIter << std::endl;
		//BOOST_CHECK( vertex ==  *simplexIter);
		//simplexIter--;
	  std::cout << "test_simplex_tree_contains - vertex=" << vertex << "||" << simplex.first.at(simplexIndex) << std::endl;
	  BOOST_CHECK(vertex ==  simplex.first.at(simplexIndex));
	  BOOST_CHECK(simplexIndex >= 0);
	  simplexIndex--;
	}
}

void test_simplex_tree_insert_returns_true(const typePairSimplexBool& returnValue)
{
	BOOST_CHECK(returnValue.second == true);
	typeST::Simplex_handle shReturned = returnValue.first; // Simplex_handle = boost::container::flat_map< Vertex_handle, Node >::iterator
	BOOST_CHECK(shReturned != typeST::Simplex_handle(NULL));
}

// Global variables
Filtration_value max_fil = 0.0;
int dim_max = -1;

void set_and_test_simplex_tree_dim_fil(typeST& simplexTree, int vectorSize, const Filtration_value& fil)
{
	if (vectorSize > dim_max+1)
	{
		dim_max = vectorSize-1;
		simplexTree.set_dimension(dim_max);
		std::cout << "   set_and_test_simplex_tree_dim_fil - dim_max=" << dim_max << std::endl;
	}
	if (fil > max_fil)
	{
		max_fil = fil;
		simplexTree.set_filtration(max_fil);
		std::cout << "   set_and_test_simplex_tree_dim_fil - max_fil=" << max_fil << std::endl;
	}
	unsigned int nb_simplices = simplexTree.num_simplices() + 1;
	simplexTree.set_num_simplices(nb_simplices);

	BOOST_CHECK( simplexTree.dimension() == dim_max );
	BOOST_CHECK( AreAlmostTheSame(simplexTree.filtration(), max_fil) );
	BOOST_CHECK( simplexTree.num_simplices() == nb_simplices );
}

BOOST_AUTO_TEST_CASE( simplex_tree_insertion )
{
	const Filtration_value FIRST_FILTRATION_VALUE = 0.1;
	const Filtration_value SECOND_FILTRATION_VALUE = 0.2;
	const Filtration_value THIRD_FILTRATION_VALUE = 0.3;
	const Filtration_value FOURTH_FILTRATION_VALUE = 0.4;
	Vertex_handle FIRST_VERTEX_HANDLE =  (Vertex_handle)0;
	Vertex_handle SECOND_VERTEX_HANDLE = (Vertex_handle)1;
	Vertex_handle THIRD_VERTEX_HANDLE = (Vertex_handle)2;
	Vertex_handle FOURTH_VERTEX_HANDLE = (Vertex_handle)3;

	// TEST OF INSERTION
	std::cout << "********************************************************************" << std::endl;
	std::cout << "TEST OF INSERTION" << std::endl;
	typeST st;

	// ++ FIRST
	std::cout << "   - INSERT 0" << std::endl;
	typeVectorVertex firstSimplexVector;
	firstSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	BOOST_CHECK( firstSimplexVector.size() == 1 );
	typeSimplex firstSimplex = std::make_pair(firstSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
	typePairSimplexBool returnValue =
			st.insert ( firstSimplex.first, firstSimplex.second );

	test_simplex_tree_insert_returns_true(returnValue);
	set_and_test_simplex_tree_dim_fil(st, firstSimplexVector.size(), firstSimplex.second);
	BOOST_CHECK( st.num_vertices() == (size_t)1 );

	// ++ SECOND
	std::cout << "   - INSERT 1" << std::endl;
	typeVectorVertex secondSimplexVector;
	secondSimplexVector.push_back(SECOND_VERTEX_HANDLE);
	BOOST_CHECK( secondSimplexVector.size() == 1 );
	typeSimplex secondSimplex = std::make_pair(secondSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
	returnValue =
			st.insert ( secondSimplex.first, secondSimplex.second );

	test_simplex_tree_insert_returns_true(returnValue);
	set_and_test_simplex_tree_dim_fil(st, secondSimplexVector.size(), secondSimplex.second);
	BOOST_CHECK( st.num_vertices() == (size_t)2 );

	// ++ THIRD
	std::cout << "   - INSERT (0,1)" << std::endl;
	typeVectorVertex thirdSimplexVector;
	thirdSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	thirdSimplexVector.push_back(SECOND_VERTEX_HANDLE);
	BOOST_CHECK( thirdSimplexVector.size() == 2 );
	typeSimplex thirdSimplex = std::make_pair(thirdSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
	returnValue =
			st.insert ( thirdSimplex.first, thirdSimplex.second );

	test_simplex_tree_insert_returns_true(returnValue);
	set_and_test_simplex_tree_dim_fil(st, thirdSimplexVector.size(), thirdSimplex.second);
	BOOST_CHECK( st.num_vertices() == (size_t)2 ); // Not incremented !!

	// ++ FOURTH
	std::cout << "   - INSERT 2" << std::endl;
	typeVectorVertex fourthSimplexVector;
	fourthSimplexVector.push_back(THIRD_VERTEX_HANDLE);
	BOOST_CHECK( fourthSimplexVector.size() == 1 );
	typeSimplex fourthSimplex = std::make_pair(fourthSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
	returnValue =
			st.insert ( fourthSimplex.first, fourthSimplex.second );

	test_simplex_tree_insert_returns_true(returnValue);
	set_and_test_simplex_tree_dim_fil(st, fourthSimplexVector.size(), fourthSimplex.second);
	BOOST_CHECK( st.num_vertices() == (size_t)3 );

	// ++ FIFTH
	std::cout << "   - INSERT (2,0)" << std::endl;
	typeVectorVertex fifthSimplexVector;
	fifthSimplexVector.push_back(THIRD_VERTEX_HANDLE);
	fifthSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	BOOST_CHECK( fifthSimplexVector.size() == 2 );
	typeSimplex fifthSimplex = std::make_pair(fifthSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
	returnValue =
			st.insert ( fifthSimplex.first, fifthSimplex.second );

	test_simplex_tree_insert_returns_true(returnValue);
	set_and_test_simplex_tree_dim_fil(st, fifthSimplexVector.size(), fifthSimplex.second);
	BOOST_CHECK( st.num_vertices() == (size_t)3 ); // Not incremented !!

	// ++ SIXTH
	std::cout << "   - INSERT (2,1)" << std::endl;
	typeVectorVertex sixthSimplexVector;
	sixthSimplexVector.push_back(THIRD_VERTEX_HANDLE);
	sixthSimplexVector.push_back(SECOND_VERTEX_HANDLE);
	BOOST_CHECK( sixthSimplexVector.size() == 2 );
	typeSimplex sixthSimplex = std::make_pair(sixthSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
	returnValue =
			st.insert ( sixthSimplex.first, sixthSimplex.second );

	test_simplex_tree_insert_returns_true(returnValue);
	set_and_test_simplex_tree_dim_fil(st, sixthSimplexVector.size(), sixthSimplex.second);
	BOOST_CHECK( st.num_vertices() == (size_t)3 ); // Not incremented !!

	// ++ SEVENTH
	std::cout << "   - INSERT (2,1,0)" << std::endl;
	typeVectorVertex seventhSimplexVector;
	seventhSimplexVector.push_back(THIRD_VERTEX_HANDLE);
	seventhSimplexVector.push_back(SECOND_VERTEX_HANDLE);
	seventhSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	BOOST_CHECK( seventhSimplexVector.size() == 3 );
	typeSimplex seventhSimplex = std::make_pair(seventhSimplexVector, Filtration_value(THIRD_FILTRATION_VALUE));
	returnValue =
			st.insert ( seventhSimplex.first, seventhSimplex.second );

	test_simplex_tree_insert_returns_true(returnValue);
	set_and_test_simplex_tree_dim_fil(st, seventhSimplexVector.size(), seventhSimplex.second);
	BOOST_CHECK( st.num_vertices() == (size_t)3 ); // Not incremented !!

	// ++ EIGHTH
	std::cout << "   - INSERT 3" << std::endl;
	typeVectorVertex eighthSimplexVector;
	eighthSimplexVector.push_back(FOURTH_VERTEX_HANDLE);
	BOOST_CHECK( eighthSimplexVector.size() == 1 );
	typeSimplex eighthSimplex = std::make_pair(eighthSimplexVector, Filtration_value(FIRST_FILTRATION_VALUE));
	returnValue =
			st.insert ( eighthSimplex.first, eighthSimplex.second );

	test_simplex_tree_insert_returns_true(returnValue);
	set_and_test_simplex_tree_dim_fil(st, eighthSimplexVector.size(), eighthSimplex.second);
	BOOST_CHECK( st.num_vertices() == (size_t)4 );

	// ++ NINETH
	std::cout << "   - INSERT (3,0)" << std::endl;
	typeVectorVertex ninethSimplexVector;
	ninethSimplexVector.push_back(FOURTH_VERTEX_HANDLE);
	ninethSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	BOOST_CHECK( ninethSimplexVector.size() == 2 );
	typeSimplex ninethSimplex = std::make_pair(ninethSimplexVector, Filtration_value(SECOND_FILTRATION_VALUE));
	returnValue =
			st.insert ( ninethSimplex.first, ninethSimplex.second );

	test_simplex_tree_insert_returns_true(returnValue);
	set_and_test_simplex_tree_dim_fil(st, ninethSimplexVector.size(), ninethSimplex.second);
	BOOST_CHECK( st.num_vertices() == (size_t)4 ); // Not incremented !!

	// ++ TENTH
	std::cout << "   - INSERT 0 (already inserted)" << std::endl;
	typeVectorVertex tenthSimplexVector;
	tenthSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	BOOST_CHECK( tenthSimplexVector.size() == 1 );
	typeSimplex tenthSimplex = std::make_pair(tenthSimplexVector, Filtration_value(FOURTH_FILTRATION_VALUE)); // With a different filtration value
	returnValue =
			st.insert ( tenthSimplex.first, tenthSimplex.second );

	BOOST_CHECK(returnValue.second == false);
	typeST::Simplex_handle shReturned = returnValue.first; // Simplex_handle = boost::container::flat_map< Vertex_handle, Node >::iterator
	BOOST_CHECK(shReturned == typeST::Simplex_handle(NULL));
	BOOST_CHECK( st.num_vertices() == (size_t)4 ); // Not incremented !!
	BOOST_CHECK( st.dimension() == dim_max );
	BOOST_CHECK( AreAlmostTheSame(st.filtration(), max_fil) );

	// ++ ELEVENTH
	std::cout << "   - INSERT (2,1,0) (already inserted)" << std::endl;
	typeVectorVertex eleventhSimplexVector;
	eleventhSimplexVector.push_back(THIRD_VERTEX_HANDLE);
	eleventhSimplexVector.push_back(SECOND_VERTEX_HANDLE);
	eleventhSimplexVector.push_back(FIRST_VERTEX_HANDLE);
	BOOST_CHECK( eleventhSimplexVector.size() == 3 );
	typeSimplex eleventhSimplex = std::make_pair(eleventhSimplexVector, Filtration_value(FOURTH_FILTRATION_VALUE));
	returnValue =
			st.insert ( eleventhSimplex.first, eleventhSimplex.second );

	BOOST_CHECK(returnValue.second == false);
	shReturned = returnValue.first; // Simplex_handle = boost::container::flat_map< Vertex_handle, Node >::iterator
	BOOST_CHECK(shReturned == typeST::Simplex_handle(NULL));
	BOOST_CHECK( st.num_vertices() == (size_t)4 ); // Not incremented !!
	BOOST_CHECK( st.dimension() == dim_max );
	BOOST_CHECK( AreAlmostTheSame(st.filtration(), max_fil) );

  /* Inserted simplex:        */
  /*    1                     */
  /*    o                     */
  /*   /X\                    */
  /*  o---o---o               */
  /*  2   0   3               */

	//   [0.1] 0
	//   [0.1] 1
	//   [0.1] 2
	//   [0.1] 3
	//   [0.2] 1 0
	//   [0.2] 2 0
	//   [0.2] 2 1
	//   [0.2] 3 0
	//   [0.3] 2 1 0

	//  !! Be careful, simplex are sorted by filtration value on insertion !!
	std::cout << "simplex_tree_insertion - first - 0" << std::endl;
	test_simplex_tree_contains(st, firstSimplex, 0); // (0) -> 0
	std::cout << "simplex_tree_insertion - second - 1" << std::endl;
	test_simplex_tree_contains(st, secondSimplex, 1); // (1) -> 1
	std::cout << "simplex_tree_insertion - third - 4" << std::endl;
	test_simplex_tree_contains(st, thirdSimplex, 4); // (0,1) -> 4
	std::cout << "simplex_tree_insertion - fourth - 2" << std::endl;
	test_simplex_tree_contains(st, fourthSimplex, 2); // (2) -> 2
	std::cout << "simplex_tree_insertion - fifth - 5" << std::endl;
	test_simplex_tree_contains(st, fifthSimplex, 5); // (2,0) -> 5
	std::cout << "simplex_tree_insertion - sixth - 6" << std::endl;
	test_simplex_tree_contains(st, sixthSimplex, 6); //(2,1) -> 6
	std::cout << "simplex_tree_insertion - seventh - 8" << std::endl;
	test_simplex_tree_contains(st, seventhSimplex, 8); // (2,1,0) -> 8
	std::cout << "simplex_tree_insertion - eighth - 3" << std::endl;
	test_simplex_tree_contains(st, eighthSimplex, 3); // (3) -> 3
	std::cout << "simplex_tree_insertion - nineth - 7" << std::endl;
	test_simplex_tree_contains(st, ninethSimplex, 7); // (3,0) -> 7

	// Display the Simplex_tree - Can not be done in the middle of 2 inserts
	std::cout << "The complex contains " << st.num_simplices() << " simplices" << std::endl;
	std::cout << "   - dimension " << st.dimension() << "   - filtration " << st.filtration() << std::endl;
	std::cout << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
	for( auto f_simplex : st.filtration_simplex_range() )
	{
		std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
		for( auto vertex : st.simplex_vertex_range(f_simplex) )
		{
			std::cout << (int)vertex << " ";
		}
		std::cout << std::endl;
	}

}
