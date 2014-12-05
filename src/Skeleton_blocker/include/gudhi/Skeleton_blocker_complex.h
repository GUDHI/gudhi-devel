/*
 * Skeleton_blocker_complex.h
 *
 *  Created on: Jan, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_SKELETON_BLOCKER_COMPLEX_H
#define GUDHI_SKELETON_BLOCKER_COMPLEX_H


#include <map>
#include <list>
#include <set>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <memory>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/adaptor/map.hpp>

#include "gudhi/Skeleton_blocker/iterators/Skeleton_blockers_iterators.h"
#include "gudhi/Skeleton_blocker_link_complex.h"
#include "gudhi/Skeleton_blocker/Skeleton_blocker_link_superior.h"
#include "gudhi/Skeleton_blocker/Skeleton_blocker_sub_complex.h"
#include "gudhi/Skeleton_blocker/Skeleton_blocker_simplex.h"

#include "gudhi/Skeleton_blocker/Skeleton_blocker_complex_visitor.h"
#include "gudhi/Skeleton_blocker/internal/Top_faces.h"
#include "gudhi/Utils.h"


namespace Gudhi {

namespace skbl {


/**
 *@class Skeleton_blocker_complex
 *@brief Abstract Simplicial Complex represented with a skeleton/blockers pair.
 *
 * A simplicial complex is completely defined by :
 * - the graph of its 1-skeleton;
 * - its set of blockers.
 *
 * The graph is a boost graph templated with SkeletonBlockerDS::Graph_vertex and SkeletonBlockerDS::Graph_edge.
 *
 * One can accesses to vertices via SkeletonBlockerDS::Vertex_handle, to edges via Skeleton_blocker_complex::Edge_handle and
 * simplices via Skeleton_blocker_simplex<Vertex_handle>.
 *
 * The SkeletonBlockerDS::Root_vertex_handle serves in the case of a subcomplex (see class Skeleton_blocker_sub_complex)
 * to access to the address of one vertex in the parent complex.
 *
 * @todo TODO Simplex_handle are not classic handle
 *
 */
template<class SkeletonBlockerDS>
class Skeleton_blocker_complex
{
	template<class ComplexType> friend class Complex_vertex_iterator;
	template<class ComplexType> friend class Complex_neighbors_vertices_iterator;
	template<class ComplexType> friend class Complex_edge_iterator;
	template<class ComplexType> friend class Complex_edge_around_vertex_iterator;


	template<class ComplexType> friend class Skeleton_blocker_link_complex;
	template<class ComplexType> friend class Skeleton_blocker_link_superior;
	template<class ComplexType> friend class Skeleton_blocker_sub_complex;

public:


	typedef typename SkeletonBlockerDS::Graph_vertex Graph_vertex;
	typedef typename SkeletonBlockerDS::Graph_edge Graph_edge;

	typedef typename SkeletonBlockerDS::Root_vertex_handle Root_vertex_handle;
	typedef typename SkeletonBlockerDS::Vertex_handle Vertex_handle;
	typedef typename Root_vertex_handle::boost_vertex_handle boost_vertex_handle;


	typedef Skeleton_blocker_simplex<Vertex_handle> Simplex_handle;
	typedef Skeleton_blocker_simplex<Root_vertex_handle> Root_simplex_handle;


	typedef Simplex_handle* Blocker_handle;



	typedef typename Root_simplex_handle::Simplex_vertex_const_iterator Root_simplex_iterator;
	typedef typename Simplex_handle::Simplex_vertex_const_iterator Simplex_handle_iterator;


protected:

	typedef typename boost::adjacency_list
			< boost::setS, //edges
			boost::vecS, // vertices
			boost::undirectedS,
			Graph_vertex,
			Graph_edge
			> Graph;
	//todo/remark : edges are not sorted, it heavily penalizes computation for SuperiorLink
	// (eg Link with greater vertices)
	// that burdens simplex iteration / complex initialization via list of simplices.
	// to avoid that, one should modify the graph by storing two lists of adjacency for every
	// vertex, the one with superior and the one with lower vertices, that way there is
	// no more extra cost for computation of SuperiorLink
	typedef typename boost::graph_traits<Graph>::vertex_iterator boost_vertex_iterator;
	typedef typename boost::graph_traits<Graph>::edge_iterator boost_edge_iterator;

protected:
	typedef typename boost::graph_traits<Graph>::adjacency_iterator boost_adjacency_iterator;

public:
	/**
	 * Handle to an edge of the complex.
	 */
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge_handle;



	typedef std::multimap<Vertex_handle,Simplex_handle *> BlockerMap;
	typedef typename std::multimap<Vertex_handle,Simplex_handle *>::value_type BlockerPair;
	typedef typename std::multimap<Vertex_handle,Simplex_handle *>::iterator BlockerMapIterator;
	typedef typename std::multimap<Vertex_handle,Simplex_handle *>::const_iterator BlockerMapConstIterator;

protected:
	int num_vertices_;
	int num_blockers_;

	typedef Skeleton_blocker_complex_visitor<Vertex_handle> Visitor;
	//	typedef Visitor* Visitor_ptr;
	Visitor* visitor;

	/**
	 * @details If 'x' is a Vertex_handle of a vertex in the complex then degree[x] = d is its degree.
	 *
	 * This quantity is updated when adding/removing edge.
	 *
	 * This is useful because the operation
	 * list.size() is done in linear time.
	 */
	std::vector<boost_vertex_handle> degree_;
	Graph skeleton; /** 1-skeleton of the simplicial complex. */


	//todo remove!!!
public:

	/** Each vertex can access to the blockers passing through it. */
	BlockerMap blocker_map_;




public:




	/////////////////////////////////////////////////////////////////////////////
	/** @name Constructors / Destructors / Initialization
	 */
	//@{
	Skeleton_blocker_complex(int num_vertices_ = 0,Visitor* visitor_=NULL):visitor(visitor_){
		clear();
		for (int i=0; i<num_vertices_; ++i){
			add_vertex();
		}
	}

private:

	/**
	 * 	this nested class is used for the next constructor that takes as an input a list of simplices.
	 * 	It stores a vector where the ith cell contains a set of i-simplices
	 *
	 */
	class Simplices_sets_from_list{
	private:
		typedef  std::set< Simplex_handle> Container_simplices;
		typedef typename Container_simplices::iterator Simplices_iterator;
		int dimension_;
		std::vector<Container_simplices > simplices_;

	public:
		Simplices_sets_from_list(std::list<Simplex_handle>& simplices):
			dimension_(-1){
			assert(!simplices.empty());

			for(auto simplex = simplices.begin() ; simplex != simplices.end(); ++simplex ){
				dimension_ = std::max(dimension_,(int)simplex->dimension());
			}
			simplices_ = std::vector<Container_simplices >(dimension_+1);

			// compute k-simplices
			for(auto simplex = simplices.begin() ; simplex != simplices.end(); ++simplex ){
				simplices_[simplex->dimension()].insert(*simplex);
			}
		}

		Simplices_iterator begin(int k){
			assert(0<= k && k<= dimension_);
			return simplices_[k].begin();
		}

		Simplices_iterator end(int k){
			assert(0<= k && k<= dimension_);
			return simplices_[k].end();
		}


		Container_simplices& simplices(int k){
			return simplices_[k];
		}

		int dimension(){
			return dimension_;
		}

		bool contains(const Simplex_handle& simplex) const{
			if(simplex.dimension()>dimension_)
				return false;
			else
				return simplices_[simplex.dimension()].find(simplex)!= simplices_[simplex.dimension()].end();
		}

		void print(){
			for(int i = 0; i < dimension_; ++i){
				std::cout << i<<"-simplices"<<std::endl;
				auto l = simplices_[i];
				for(auto s : l){
					std::cout << s<<std::endl;
				}
			}
		}
	};


	void compute_next_expand(
			Simplices_sets_from_list& simplices,
			int dim,
			std::list<Simplex_handle>& next_expand)
	{
		next_expand.clear();

		//	    high_resolution_clock::time_point tbeginfor = high_resolution_clock::now();
		//	    auto durationlooplink = std::chrono::duration_cast<std::chrono::microseconds>( tbeginfor - tbeginfor ).count();

		for(auto sigma = simplices.begin(dim); sigma != simplices.end(dim); ++sigma){
			//		    high_resolution_clock::time_point tbeg = high_resolution_clock::now();
			Simplex_handle t(*sigma);
			Skeleton_blocker_link_superior<Skeleton_blocker_complex> link(*this,t);
			// xxx all time here, most likely because accessing superior edges needs passing through lower one
			// currently

			//			durationlooplink += std::chrono::duration_cast<std::chrono::microseconds>( high_resolution_clock::now() - tbeg ).count();

			for(auto v : link.vertex_range()){
				Vertex_handle v_in_complex(*this->get_address(  link.get_id(v)) );
				t.add_vertex(v_in_complex);
				next_expand.push_back(t);
				t.remove_vertex(v_in_complex);
			}
		}
		//	    high_resolution_clock::time_point t2 = high_resolution_clock::now();
		//	    auto durationlooptotal = std::chrono::duration_cast<std::chrono::microseconds>( t2 - tbeginfor ).count();
		//	    DBGVALUE(durationlooptotal);
		//	    DBGVALUE(durationlooplink);
	}



public:
	/**
	 * @brief Constructor with a list of simplices
	 * @details The list of simplices must be the list
	 * of simplices of a simplicial complex.
	 *todo rewrite as iterators
	 */
	Skeleton_blocker_complex(std::list<Simplex_handle>& simplices,Visitor* visitor_=NULL):
		num_vertices_(0),num_blockers_(0),
		visitor(visitor_){
		Simplices_sets_from_list set_simplices(simplices);

		int dim = set_simplices.dimension();


		// add 1-skeleton to the complex
		for(auto v_it = set_simplices.begin(0); v_it != set_simplices.end(0); ++v_it)
			add_vertex();

		for(auto e_it = set_simplices.begin(1); e_it != set_simplices.end(1); ++e_it){
			Vertex_handle a = e_it->first_vertex();
			Vertex_handle b = e_it->last_vertex();
			assert(contains_vertex(a) && contains_vertex(b));
			add_edge(a,b);
		}

		// then add blockers
		for(int current_dim = 1 ; current_dim <=dim ; ++current_dim){
			std::list<Simplex_handle> expansion_simplices;
			compute_next_expand(set_simplices,current_dim,expansion_simplices);

			for(const auto &simplex : expansion_simplices) {
				if(!set_simplices.contains(simplex)){
					add_blocker(simplex);
				}
			}
		}
	}



	// We cannot use the default copy constructor since we need
	// to make a copy of each of the blockers
	Skeleton_blocker_complex(const Skeleton_blocker_complex& copy){
		visitor = NULL;
		degree_ = copy.degree_;
		skeleton = Graph(copy.skeleton);
		num_vertices_ = copy.num_vertices_;

		num_blockers_ = 0;
		// we copy the blockers
		for (auto blocker : copy.const_blocker_range()){
			add_blocker(*blocker);
		}
	}

	Skeleton_blocker_complex& operator=(const Skeleton_blocker_complex& copy){
		clear();
		visitor = NULL;
		degree_ = copy.degree_;
		skeleton = Graph(copy.skeleton);
		num_vertices_ = copy.num_vertices_;

		num_blockers_ = 0;
		// we copy the blockers
		for (auto blocker : copy.const_blocker_range())
			add_blocker(*blocker);
		return *this;
	}


	/**
	 * The destructor delete all blockers allocated.
	 */
	virtual ~Skeleton_blocker_complex(){
		clear();
	}

	/**
	 * Clears the simplicial complex. After a call to this function,
	 * blockers are destroyed. The 1-skeleton and the set of blockers
	 * are both empty.
	 */
	virtual void clear(){
		// xxx for now the responsabilty of freeing the visitor is for
		// the user
		visitor = NULL;

		degree_.clear();
		num_vertices_ =0;

		remove_blockers();

		skeleton.clear();
	}

	void set_visitor(Visitor* other_visitor){
		visitor = other_visitor;
	}

	//@}




	/////////////////////////////////////////////////////////////////////////////
	/** @name Vertices operations
	 */
	//@{

public:

	/**
	 * @brief Return a local Vertex_handle of a vertex given a global one.
	 * @remark Assume that the vertex is present in the complex.
	 */
	Vertex_handle operator[](Root_vertex_handle global) const{
		auto local(get_address(global));
		assert(local);
		return *local;
	}

	Graph_vertex& operator[](Vertex_handle address){
		assert(0<=address.vertex && address.vertex< boost::num_vertices(skeleton));
		return skeleton[address.vertex];
	}

	const Graph_vertex& operator[](Vertex_handle address) const{
		assert(0<=address.vertex && address.vertex< boost::num_vertices(skeleton));
		return skeleton[address.vertex];
	}

	/**
	 * @brief Adds a vertex to the simplicial complex and returns its Vertex_handle.
	 */
	Vertex_handle add_vertex(){
		Vertex_handle address(boost::add_vertex(skeleton));
		num_vertices_++;
		(*this)[address].activate();
		// safe since we now that we are in the root complex and the field 'address' and 'id'
		// are identical for every vertices
		(*this)[address].set_id(Root_vertex_handle(address.vertex));
		degree_.push_back(0);
		if (visitor) visitor->on_add_vertex(address);
		return address;
	}


	/**
	 * @brief Remove a vertex from the simplicial complex
	 * @remark In fact, it just deactivates the vertex.
	 */
	void remove_vertex(Vertex_handle address){
		assert(contains_vertex(address));
		// We remove b
		boost::clear_vertex(address.vertex,skeleton);
		(*this)[address].deactivate();
		num_vertices_--;
		degree_[address.vertex]=-1;
		if (visitor) visitor->on_remove_vertex(address);
	}

	/**
	 * @return true iff the simplicial complex contains the vertex u
	 */
	bool contains_vertex(Vertex_handle u) const{
		if (u.vertex<0 || u.vertex>=boost::num_vertices(skeleton)) return false;
		return (*this)[u].is_active();
	}

	/**
	 * @return true iff the simplicial complex contains the vertex u
	 */
	bool contains_vertex(Root_vertex_handle u) const{
		boost::optional<Vertex_handle> address = get_address(u);
		return address &&  (*this)[*address].is_active();
	}


	/**
	 * @return true iff the simplicial complex contains all vertices
	 * of simplex sigma
	 */
	bool contains_vertices(const Simplex_handle & sigma) const{
		for (auto vertex : sigma)
			if(!contains_vertex(vertex)) return false;
		return true;
	}

	/**
	 * Given an Id return the address of the vertex having this Id in the complex.
	 * For a simplicial complex, the address is the id but it may not be the case for a SubComplex.
	 *
	 */
	virtual boost::optional<Vertex_handle> get_address(Root_vertex_handle id) const{
		boost::optional<Vertex_handle> res;
		if ( id.vertex< boost::num_vertices(skeleton) ) res = Vertex_handle(id.vertex);//xxx
		return res;
	}



	/**
	 * return the id of a vertex of adress local present in the graph
	 */
	Root_vertex_handle get_id(Vertex_handle local) const{
		assert(0<=local.vertex && local.vertex< boost::num_vertices(skeleton));
		return (*this)[local].get_id();
	}


	/**
	 * If the current complex is a sub (or sup) complex of 'other', it converts
	 * the address of a vertex v expressed in 'other' to the address of the vertex
	 * v in the current one.
	 * @remark this methods uses Root_vertex_handle to identify the vertex
	 */
	Vertex_handle convert_handle_from_another_complex(
			const Skeleton_blocker_complex& other,Vertex_handle vh_in_other) const{
		auto vh_in_current_complex = get_address(other.get_id(vh_in_other));
		assert(vh_in_current_complex);
		return *vh_in_current_complex;
	}


	int degree(Vertex_handle local) const{
		assert(0<=local.vertex && local.vertex< boost::num_vertices(skeleton));
		return degree_[local.vertex];
	}

	//@}

	/////////////////////////////////////////////////////////////////////////////
	/** @name Edges operations
	 */
	//@{

public:

	/**
	 * @brief return an edge handle if the two vertices forms
	 * an edge in the complex
	 */
	boost::optional<Edge_handle> operator[](const std::pair<Vertex_handle,Vertex_handle>& ab) const{
		boost::optional<Edge_handle> res;
		std::pair<Edge_handle,bool> edge_pair(boost::edge(ab.first.vertex,ab.second.vertex,skeleton));
		if (edge_pair.second)
			res = edge_pair.first;
		return res;
	}

	Graph_edge& operator[](Edge_handle edge_handle){
		return skeleton[edge_handle];
	}

	const Graph_edge& operator[](Edge_handle edge_handle) const{
		return skeleton[edge_handle];
	}

	Vertex_handle first_vertex(Edge_handle edge_handle) const{
		return source(edge_handle,skeleton);
	}

	Vertex_handle second_vertex(Edge_handle edge_handle) const{
		return target(edge_handle,skeleton);
	}

	/**
	 * @brief returns the simplex made with the two vertices of the edge
	 */
	Simplex_handle get_vertices(Edge_handle edge_handle) const{
		auto edge((*this)[edge_handle]);
		return Simplex_handle((*this)[edge.first()],(*this)[edge.second()]);
	}

	/**
	 * @brief Adds an edge between vertices a and b
	 */
	Edge_handle add_edge(Vertex_handle a, Vertex_handle b){
		assert(contains_vertex(a) && contains_vertex(b));
		assert(a!=b);

		auto edge_handle((*this)[std::make_pair(a,b)]);
		//		std::pair<Edge_handle,bool> pair_descr_bool = (*this)[std::make_pair(a,b)];
		//		Edge_handle edge_descr;
		//		bool edge_present = pair_descr_bool.second;
		if (!edge_handle)
		{
			edge_handle = boost::add_edge(a.vertex,b.vertex,skeleton).first;
			(*this)[*edge_handle].setId(get_id(a),get_id(b));
			degree_[a.vertex]++;
			degree_[b.vertex]++;
			if (visitor) visitor->on_add_edge(a,b);
		}
		return *edge_handle;
	}

	/**
	 * @brief Adds all edges of simplex sigma to the simplicial complex.
	 */
	void add_edges(const Simplex_handle & sigma){
		Simplex_handle_iterator i, j;
		for (i = sigma.begin() ; i != sigma.end() ; ++i)
			for (j = i, j++ ; j != sigma.end() ; ++j)
				add_edge(*i,*j);
	}

	/**
	 * @brief Removes edge ab from the simplicial complex.
	 */
	virtual Edge_handle remove_edge(Vertex_handle a, Vertex_handle b){
		bool found;
		Edge_handle edge;
		tie(edge,found) = boost::edge(a.vertex,b.vertex,skeleton);
		if (found)
		{
			if (visitor) visitor->on_remove_edge(a,b);
			//		if (heapCollapse.Contains(edge)) heapCollapse.Delete(edge);
			boost::remove_edge(a.vertex,b.vertex,skeleton);
			degree_[a.vertex]--;
			degree_[b.vertex]--;
		}
		return edge;
	}


	/**
	 * @brief Removes edge from the simplicial complex.
	 */
	void remove_edge(Edge_handle edge){
		assert(contains_vertex(first_vertex(edge)));
		assert(contains_vertex(second_vertex(edge)));
		remove_edge(first_vertex(edge),second_vertex(edge));
	}


	/**
	 * @brief The complex is reduced to its set of vertices.
	 * All the edges and blockers are removed.
	 */
	void keep_only_vertices(){
		remove_blockers();

		for(auto u : vertex_range()){
			while (this->degree(u)> 0)
			{
				Vertex_handle v(*(adjacent_vertices(u.vertex, this->skeleton).first));
				this->remove_edge(u,v);
			}
		}
	}


	/**
	 * @return true iff the simplicial complex contains an edge between
	 * vertices a and b
	 */
	bool contains_edge(Vertex_handle a, Vertex_handle b) const{
		//if (a.vertex<0 || b.vertex <0) return false;
		return boost::edge(a.vertex,b.vertex,skeleton).second;
	}


	/**
	 * @return true iff the simplicial complex contains all vertices
	 * and all edges of simplex sigma
	 */
	bool contains_edges(const Simplex_handle & sigma) const{
		for (auto i = sigma.begin() ; i != sigma.end() ; ++i){
			if(!contains_vertex(*i)) return false;
			for (auto j=i; ++j != sigma.end() ; ){
				if (!contains_edge(*i,*j))
					return false;
			}
		}
		return true;
	}
	//@}



	/////////////////////////////////////////////////////////////////////////////
	/** @name Blockers operations
	 */
	//@{
	/**
	 * Adds the 2-blocker abc
	 */
	void add_blocker(Vertex_handle a, Vertex_handle b, Vertex_handle c){
		add_blocker(Simplex_handle(a,b,c));
	}

	/**
	 * Adds the 3-blocker abcd
	 */
	void add_blocker(Vertex_handle a, Vertex_handle b, Vertex_handle c, Vertex_handle d){
		add_blocker(Simplex_handle(a,b,c,d));
	}

	/**
	 * Adds the simplex blocker_pt to the set of blockers and
	 * returns a Blocker_handle toward it if was not present before.
	 */
	Blocker_handle add_blocker(const Simplex_handle& blocker){
		if (contains_blocker(blocker))
		{
			//std::cerr << "ATTEMPT TO ADD A BLOCKER ALREADY THERE ---> BLOCKER IGNORED" << endl;
			return 0;
		}
		else{
			if (visitor) visitor->on_add_blocker(blocker);
			Blocker_handle blocker_pt = new Simplex_handle(blocker);
			num_blockers_++;
			auto vertex = blocker_pt->begin();
			while(vertex != blocker_pt->end())
			{
				blocker_map_.insert(BlockerPair(*vertex,blocker_pt));
				++vertex;
			}
			return blocker_pt;
		}
	}

protected:
	/**
	 * Adds the simplex s to the set of blockers
	 */
	void add_blocker(Blocker_handle blocker){
		if (contains_blocker(*blocker))
		{
			//std::cerr << "ATTEMPT TO ADD A BLOCKER ALREADY THERE ---> BLOCKER IGNORED" << endl;
			return;
		}
		else{
			if (visitor) visitor->on_add_blocker(*blocker);
			num_blockers_++;
			auto vertex = blocker->begin();
			while(vertex != blocker->end())
			{
				blocker_map_.insert(BlockerPair(*vertex,blocker));
				++vertex;
			}
		}
	}


protected:
	/**
	 * Removes sigma from the blocker map of vertex v
	 */
	void remove_blocker(const Blocker_handle sigma, Vertex_handle v){
		Complex_blocker_around_vertex_iterator blocker;
		for (blocker = blocker_range(v).begin();
				blocker != blocker_range(v).end();
				++blocker
		){
			if (*blocker == sigma) break;
		}
		if (*blocker != sigma){
			std::cerr << "bug ((*blocker).second == sigma) ie try to remove a blocker not present\n";
			assert(false);
		}
		else{
			blocker_map_.erase(blocker.current_position());
		}
	}

public:
	/**
	 * Removes the simplex sigma from the set of blockers.
	 * sigma has to belongs to the set of blockers
	 */
	void remove_blocker(const Blocker_handle sigma){
		for (auto vertex : *sigma){
			remove_blocker(sigma,vertex);
		}
		num_blockers_--;
	}




	/**
	 * Remove all blockers, in other words, it expand the simplicial
	 * complex to the smallest flag complex that contains it.
	 */
	void remove_blockers(){
		// Desallocate the blockers
		while (!blocker_map_.empty()){
			delete_blocker(blocker_map_.begin()->second);
		}
		num_blockers_ = 0;
		blocker_map_.clear();
	}

protected:
	/**
	 * Removes the simplex sigma from the set of blockers.
	 * sigma has to belongs to the set of blockers
	 *
	 * @remark contrarily to delete_blockers does not call the destructor
	 */
	void remove_blocker(const Simplex_handle& sigma){
		assert(contains_blocker(sigma));
		for (auto vertex : sigma)
			remove_blocker(sigma,vertex);
		num_blockers_--;
	}




public:
	/**
	 * Removes the simplex s from the set of blockers
	 * and desallocate s.
	 */
	void delete_blocker(Blocker_handle sigma){
		if (visitor) visitor->on_delete_blocker(sigma);
		remove_blocker(sigma);
		delete sigma;
	}


	/**
	 * @return true iff s is a blocker of the simplicial complex
	 */
	bool contains_blocker(const Blocker_handle s) const{
		if (s->dimension()<2)
			return false;

		Vertex_handle a = s->first_vertex();

		for (auto blocker : const_blocker_range(a)){
			if ( s == *blocker )
				return true;
		}
		return false;
	}

	/**
	 * @return true iff s is a blocker of the simplicial complex
	 */
	bool contains_blocker(const Simplex_handle & s) const{
		if (s.dimension()<2)
			return false;

		Vertex_handle a = s.first_vertex();

		for (auto blocker : const_blocker_range(a)){
			if ( s == *blocker )
				return true;
		}
		return false;
	}


private:
	/**
	 * @return true iff a blocker of the simplicial complex
	 * is a face of sigma.
	 */
	bool blocks(const Simplex_handle & sigma) const{

		for(auto blocker : const_blocker_range()){
			if ( sigma.contains(*blocker) )
				return true;
		}
		return false;
	}

	//@}




	/////////////////////////////////////////////////////////////////////////////
	/** @name Neighbourhood access
	 */
	//@{

public:
	/**
	 * @brief Adds to simplex n the neighbours of v:
	 * \f$ n \leftarrow n \cup N(v) \f$.
	 *
	 * If keep_only_superior is true then only vertices that are greater than v are added.
	 */
	virtual void add_neighbours(Vertex_handle v, Simplex_handle & n,bool keep_only_superior=false) const{
		boost_adjacency_iterator ai, ai_end;
		for (tie(ai, ai_end) = adjacent_vertices(v.vertex, skeleton); ai != ai_end; ++ai){
			if (keep_only_superior){
				if (*ai>v.vertex)
					n.add_vertex(Vertex_handle(*ai));
			}
			else
				n.add_vertex(Vertex_handle(*ai));
		}
	}

	/**
	 * @brief Add to simplex res all vertices which are
	 * neighbours of alpha: ie \f$ res \leftarrow res \cup N(alpha) \f$.
	 *
	 * If 'keep_only_superior' is true then only vertices that are greater than alpha are added.
	 *
	 * todo revoir
	 *
	 */
	virtual void add_neighbours(const Simplex_handle &alpha, Simplex_handle & res,bool keep_only_superior=false) const{
		res.clear();
		// ----------------------------
		// Compute vertices in the link
		// we compute the intersection of N(alpha_i) and store it in n
		// ----------------------------
		auto alpha_vertex = alpha.begin();
		add_neighbours(*alpha_vertex,res,keep_only_superior);
		for (alpha_vertex = (alpha.begin())++ ; alpha_vertex != alpha.end() ; ++alpha_vertex)
		{
			keep_neighbours(*alpha_vertex,res,keep_only_superior);
		}
	}

	/**
	 * @brief Eliminates from simplex n all vertices which are
	 * not neighbours of v: \f$ res \leftarrow res \cap N(v) \f$.
	 *
	 * If 'keep_only_superior' is true then only vertices that are greater than v are keeped.
	 *
	 */
	virtual void keep_neighbours(Vertex_handle v, Simplex_handle& res,bool keep_only_superior=false) const{
		Simplex_handle nv;
		add_neighbours(v,nv,keep_only_superior);
		res.intersection(nv);
	}

	/**
	 * @brief Eliminates from simplex n all vertices which are
	 * neighbours of v: \f$ res \leftarrow res \setminus N(v) \f$.
	 *
	 * If 'keep_only_superior' is true then only vertices that are greater than v are added.
	 *
	 */
	virtual void remove_neighbours(Vertex_handle v, Simplex_handle & res,bool keep_only_superior=false) const{
		Simplex_handle nv;
		add_neighbours(v,nv,keep_only_superior);
		res.difference(nv);
	}
	//@}


	/////////////////////////////////////////////////////////////////////////////
	/** @name Operations on the simplicial complex
	 */
	//@{
public:

	/**
	 * @brief Compute the local vertices of 's' in the current complex
	 * If one of them is not present in the complex then the return value is uninitialized.
	 *
	 * xxx rename get_address et place un using dans sub_complex
	 */
	boost::optional<Simplex_handle>	get_simplex_address(const Root_simplex_handle& s) const
	{
		boost::optional<Simplex_handle> res;

		Simplex_handle s_address;
		//Root_simplex_const_iterator i;
		for (auto i = s.begin() ; i != s.end() ; ++i)
		{
			boost::optional<Vertex_handle> address = get_address(*i);
			if (!address)
				return res;
			else
				s_address.add_vertex(*address);
		}
		res = s_address;
		return res;
	}

	/**
	 * @brief returns a simplex with vertices which are the id of vertices of the
	 * argument.
	 */
	Root_simplex_handle get_id(const Simplex_handle& local_simplex) const{
		Root_simplex_handle global_simplex;
		for (auto x = local_simplex.begin(); x!= local_simplex.end();++x){
			global_simplex.add_vertex(get_id(*x));

		}
		return global_simplex;

	}


	/**
	 * @brief returns true iff the simplex s belongs to the simplicial
	 * complex.
	 */
	virtual bool contains(const Simplex_handle & s) const{
		if (s.dimension() == -1 ) return false;
		else
			if (s.dimension() ==0 ){
				return contains_vertex(s.first_vertex());
			}
			else
				return ( contains_edges(s) && !blocks(s) );
	}

	/*
	 * @brief returnrs true iff the complex is empty.
	 */
	bool empty() const{
		return num_vertices()==0;
	}

	/*
	 * @brief returns the number of vertices in the complex.
	 */
	int num_vertices() const{
		//remark boost::num_vertices(skeleton) counts deactivated vertices
		return num_vertices_;
	}

	/*
	 * @brief returns the number of edges in the complex.
	 * todo in O(n), cache the value
	 */
	int num_edges() const{
		return boost::num_edges(skeleton);
	}

	/*
	 * @brief returns the number of blockers in the complex.
	 */
	int num_blockers() const{
		return num_blockers_;
	}

	/*
	 * @brief returns true iff the graph of the 1-skeleton of the complex is complete.
	 */
	bool complete() const{
		return (num_vertices()*(num_vertices()-1))/2 == num_edges();
	}

	/**
	 * @brief returns the number of connected components in the graph of the 1-skeleton.
	 */
	int num_connected_components() const{
		int num_vert_collapsed = skeleton.vertex_set().size() - num_vertices();
		std::vector<int> component(skeleton.vertex_set().size());
		return boost::connected_components(this->skeleton,&component[0]) - num_vert_collapsed;
	}


	//todo remove
	// do
	void keep_only_largest_cc(){
		std::vector<unsigned> component(skeleton.vertex_set().size());
		boost::connected_components(this->skeleton,&component[0]);
		auto maxCC = min_element(component.begin(),component.end());
		for (unsigned i = 0; i != component.size(); ++i){
			if(component[i]!=*maxCC){
				if(this->contains_vertex(Vertex_handle(i)))
					this->remove_vertex(Vertex_handle(i));
			}
		}
	}


	/**
	 * @brief %Test if the complex is a cone.
	 * @details Runs in O(n) where n is the number of vertices.
	 */
	bool is_cone() const{
		if (num_vertices()==0) return false;
		if (num_vertices()==1) return true;
		for(auto vi :  vertex_range()){
			//xxx todo faire une methode bool is_in_blocker(Vertex_handle)
			if (blocker_map_.find(vi)==blocker_map_.end()){
				// no blocker passes through the vertex, we just need to
				// check if the current vertex is linked to all others vertices of the complex
				if (degree_[vi.vertex] == num_vertices()-1)
					return true;
			}
		}
		return false;
	}

	//@}


	/////////////////////////////////////////////////////////////////////////////
	/** @name Vertex iterators
	 */
	//@{

	typedef Complex_vertex_iterator<Skeleton_blocker_complex> CVI; //todo rename

	//	/**
	//	 * @brief Range over the vertices of the simplicial complex.
	//	 * Methods .begin() and .end() return a Complex_vertex_iterator.
	//	 */
	typedef boost::iterator_range     < Complex_vertex_iterator<Skeleton_blocker_complex> > Complex_vertex_range;

	/**
	 * @brief Returns a Complex_vertex_range over all vertices of the complex
	 */
	Complex_vertex_range vertex_range() const
	{
		auto begin = Complex_vertex_iterator<Skeleton_blocker_complex>(this);
		auto end = Complex_vertex_iterator<Skeleton_blocker_complex>(this,0);
		return Complex_vertex_range(begin,end);
	}

	typedef boost::iterator_range     < Complex_neighbors_vertices_iterator<Skeleton_blocker_complex> > Complex_neighbors_vertices_range;

	/**
	 * @brief Returns a Complex_edge_range over all edges of the simplicial complex that passes trough v
	 */
	Complex_neighbors_vertices_range vertex_range(Vertex_handle v) const
	{
		auto begin = Complex_neighbors_vertices_iterator<Skeleton_blocker_complex>(this,v);
		auto end = Complex_neighbors_vertices_iterator<Skeleton_blocker_complex>(this,v,0);
		return Complex_neighbors_vertices_range(begin,end);
	}

	//@}


	/** @name Edge iterators
	 */
	//@{


	typedef boost::iterator_range  <Complex_edge_iterator<Skeleton_blocker_complex<SkeletonBlockerDS>>>
			Complex_edge_range;

	/**
	 * @brief Returns a Complex_edge_range over all edges of the simplicial complex
	 */
	Complex_edge_range edge_range() const
	{
		auto begin = Complex_edge_iterator<Skeleton_blocker_complex<SkeletonBlockerDS>>(this);
		auto end = Complex_edge_iterator<Skeleton_blocker_complex<SkeletonBlockerDS>>(this,0);
		return Complex_edge_range(begin,end);
	}



	typedef boost::iterator_range  <Complex_edge_around_vertex_iterator<Skeleton_blocker_complex<SkeletonBlockerDS>>>
			Complex_edge_around_vertex_range;
	/**
	 * @brief Returns a Complex_edge_range over all edges of the simplicial complex that passes
	 * through 'v'
	 */
	Complex_edge_around_vertex_range edge_range(Vertex_handle v) const
	{
		auto begin = Complex_edge_around_vertex_iterator<Skeleton_blocker_complex<SkeletonBlockerDS>>(this,v);
		auto end = Complex_edge_around_vertex_iterator<Skeleton_blocker_complex<SkeletonBlockerDS>>(this,v,0);
		return Complex_edge_around_vertex_range(begin,end);
	}



	//@}


	/** @name Triangles iterators
	 */
	//@{
private:
	typedef Skeleton_blocker_link_complex<Skeleton_blocker_complex<SkeletonBlockerDS> > Link;
	typedef Skeleton_blocker_link_superior<Skeleton_blocker_complex<SkeletonBlockerDS> > Superior_link;
public:
	typedef Triangle_around_vertex_iterator<Skeleton_blocker_complex,Superior_link> Superior_triangle_around_vertex_iterator;

	typedef boost::iterator_range     < Triangle_around_vertex_iterator<Skeleton_blocker_complex,Link> > Complex_triangle_around_vertex_range;

	/**
	 * @brief Range over triangles around a vertex of the simplicial complex.
	 * Methods .begin() and .end() return a Triangle_around_vertex_iterator.
	 *
	 */
	Complex_triangle_around_vertex_range triangle_range(Vertex_handle v) const
	{
		auto begin = Triangle_around_vertex_iterator<Skeleton_blocker_complex,Link>(this,v);
		auto end = Triangle_around_vertex_iterator<Skeleton_blocker_complex,Link>(this,v,0);
		return Complex_triangle_around_vertex_range(begin,end);
	}


	typedef boost::iterator_range<Triangle_iterator<Skeleton_blocker_complex> > Complex_triangle_range;

	/**
	 * @brief Range over triangles of the simplicial complex.
	 * Methods .begin() and .end() return a Triangle_around_vertex_iterator.
	 *
	 */
	Complex_triangle_range triangle_range() const
	{
		auto end = Triangle_iterator<Skeleton_blocker_complex>(this,0);
		if(empty()){
			return Complex_triangle_range(end,end);
		}
		else{
			auto begin = Triangle_iterator<Skeleton_blocker_complex>(this);
			return Complex_triangle_range(begin,end);
		}

	}


	//@}



	/** @name Simplices iterators
	 */
	//@{
	typedef Simplex_around_vertex_iterator<Skeleton_blocker_complex,Link> Complex_simplex_around_vertex_iterator;

	/**
	 * @brief Range over the simplices of the simplicial complex around a vertex.
	 * Methods .begin() and .end() return a Complex_simplex_around_vertex_iterator.
	 */
	typedef boost::iterator_range     < Complex_simplex_around_vertex_iterator > Complex_simplex_around_vertex_range;

	/**
	 * @brief Returns a Complex_simplex_around_vertex_range over all the simplices around a vertex of the complex
	 */
	Complex_simplex_around_vertex_range simplex_range(Vertex_handle v) const
	{
		assert(contains_vertex(v));
		return Complex_simplex_around_vertex_range(
				Complex_simplex_around_vertex_iterator(this,v),
				Complex_simplex_around_vertex_iterator(this,v,true)
		);
	}

	//	typedef Simplex_iterator<Skeleton_blocker_complex,Superior_link> Complex_simplex_iterator;
	typedef Simplex_iterator<Skeleton_blocker_complex> Complex_simplex_iterator;

	typedef boost::iterator_range     < Complex_simplex_iterator > Complex_simplex_range;

	Complex_simplex_range simplex_range() const
	{
		Complex_simplex_iterator end(this,true);
		if(empty()){
			return Complex_simplex_range(end,end);
		}
		else{
			Complex_simplex_iterator begin(this);
			return Complex_simplex_range(begin,end);
		}
	}


	//@}


	/** @name Blockers iterators
	 */
	//@{
private:
	/**
	 * @brief Iterator over the blockers adjacent to a vertex
	 */
	typedef Blocker_iterator_around_vertex_internal<
			typename std::multimap<Vertex_handle,Simplex_handle *>::iterator,
			Blocker_handle>
	Complex_blocker_around_vertex_iterator;

	/**
	 * @brief Iterator over (constant) blockers adjacent to a vertex
	 */
	typedef Blocker_iterator_around_vertex_internal<
			typename std::multimap<Vertex_handle,Simplex_handle *>::const_iterator,
			const Blocker_handle>
	Const_complex_blocker_around_vertex_iterator;

	typedef boost::iterator_range  <Complex_blocker_around_vertex_iterator>	Complex_blocker_around_vertex_range;
	typedef boost::iterator_range  <Const_complex_blocker_around_vertex_iterator>	Const_complex_blocker_around_vertex_range;

public:



	Complex_blocker_around_vertex_range blocker_range(Vertex_handle v)
	{
		auto begin = Complex_blocker_around_vertex_iterator(blocker_map_.lower_bound(v));
		auto end = Complex_blocker_around_vertex_iterator(blocker_map_.upper_bound(v));
		return Complex_blocker_around_vertex_range(begin,end);
	}

	/**
	 * @brief Returns a Complex_blocker_around_vertex_range over all blockers of the complex adjacent to the vertex v.
	 */
	Const_complex_blocker_around_vertex_range const_blocker_range(Vertex_handle v) const
	{
		auto begin = Const_complex_blocker_around_vertex_iterator(blocker_map_.lower_bound(v));
		auto end = Const_complex_blocker_around_vertex_iterator(blocker_map_.upper_bound(v));
		return Const_complex_blocker_around_vertex_range(begin,end);
	}




private:

	/**
	 * @brief Iterator over the blockers.
	 */
	typedef Blocker_iterator_internal<
			typename std::multimap<Vertex_handle,Simplex_handle *>::iterator,
			Blocker_handle>
	Complex_blocker_iterator;

	/**
	 * @brief Iterator over the (constant) blockers.
	 */
	typedef Blocker_iterator_internal<
			typename std::multimap<Vertex_handle,Simplex_handle *>::const_iterator,
			const Blocker_handle>
	Const_complex_blocker_iterator;

	typedef boost::iterator_range  <Complex_blocker_iterator>	Complex_blocker_range;
	typedef boost::iterator_range  <Const_complex_blocker_iterator>	Const_complex_blocker_range;


public:


	Complex_blocker_range blocker_range()
	{
		auto begin = Complex_blocker_iterator(blocker_map_.begin(), blocker_map_.end() );
		auto end = Complex_blocker_iterator(blocker_map_.end() , blocker_map_.end() );
		return Complex_blocker_range(begin,end);
	}


	Const_complex_blocker_range const_blocker_range() const
	{
		auto begin = Const_complex_blocker_iterator(blocker_map_.begin(), blocker_map_.end() );
		auto end = Const_complex_blocker_iterator(blocker_map_.end() , blocker_map_.end() );
		return Const_complex_blocker_range(begin,end);
	}




	//@}





	/////////////////////////////////////////////////////////////////////////////
	/** @name Print and IO methods
	 */
	//@{
public:
	std::string to_string() const{
		std::ostringstream stream;
		stream<<num_vertices()<<" vertices:\n"<<vertices_to_string()<<std::endl;
		stream<<num_edges()<<" edges:\n"<<edges_to_string()<<std::endl;
		stream<<num_blockers()<<" blockers:\n"<<blockers_to_string()<<std::endl;
		return stream.str();
	}

	std::string vertices_to_string() const{
		std::ostringstream stream;
		for(auto vertex : vertex_range())
			stream << "("<<(*this)[vertex].get_id()<<"),";
		stream<< std::endl;
		return stream.str();
	}

	std::string edges_to_string() const{
		std::ostringstream stream;
		for(auto edge : edge_range())
			stream << "("<< (*this)[edge].first()<<","<< (*this)[edge].second() << ")"<<" id = "<< (*this)[edge].index()<< std::endl;
		stream<< std::endl;
		return stream.str();
	}


	std::string blockers_to_string() const{
		std::ostringstream stream;
		for (auto bl:blocker_map_){
			stream << bl.first << " => " << bl.second << ":"<<*bl.second <<"\n";
		}
		return stream.str();
	}

	//@}

};



/**
 * build a simplicial complex from a collection
 * of top faces.
 */
template<typename Complex,typename SimplexHandleIterator>
unsigned make_complex_from_top_faces(Complex& complex,SimplexHandleIterator begin,SimplexHandleIterator end){
	typedef typename Complex::Simplex_handle Simplex_handle;

	int dimension = 0;
	for(auto top_face = begin; top_face != end; ++top_face)
		dimension = std::max(dimension,top_face->dimension());

	std::vector< std::set<Simplex_handle> > simplices_per_dimension(dimension+1);


	for(auto top_face = begin; top_face != end; ++top_face){
		register_faces(simplices_per_dimension,*top_face);
	}

	// compute list of simplices
	std::list<Simplex_handle> simplices;
	for(int dim = 0 ; dim <= dimension ; ++dim){
		std::copy(
				simplices_per_dimension[dim].begin(),
				simplices_per_dimension[dim].end(),
				std::back_inserter(simplices)
		);
	}
	complex = Complex(simplices);
	return simplices.size();
}



}

}  // namespace GUDHI



#endif

