 /*    This file is part of the Gudhi Library. The Gudhi library 
  *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
  *    library for computational topology.
  *
  *    Author(s):       Clément Maria
  *
  *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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

#ifndef _PERSISTENCECOMPUTATION_SIMPLEXTREE_
#define _PERSISTENCECOMPUTATION_SIMPLEXTREE_

#include <boost/tuple/tuple.hpp>
#include <boost/intrusive/set.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/intrusive/list.hpp>
#include <boost/pool/object_pool.hpp>
#include "gudhi/Persistent_cohomology/Persistent_cohomology_column.h"
#include "gudhi/Persistent_cohomology/Field_Zp.h"

namespace Gudhi{

/** \defgroup persistent_cohomology Persistent Cohomology
  *

  Computation of persistent cohomology using the algorithm of 
   \cite DBLP:journals/dcg/SilvaMV11 and \cite DBLP:journals/corr/abs-1208-5018 
   and the Compressed Annotation Matrix 
   implementation of \cite DBLP:conf/esa/BoissonnatDM13 
       
  The theory of homology consists in attaching to a topological space a sequence of 
  (homology) groups, 
  capturing global topological features 
  like connected components, holes, cavities, etc. Persistent homology studies the evolution 
  -- birth, life and death -- of 
  these features when the topological space is changing. Consequently, the theory is essentially 
  composed of three elements: 
  topological spaces, their homology groups and an evolution scheme.

  <DT>Topological Spaces:</DT>
  Topological spaces are represented by simplicial complexes.
  Let \f$V = \{1, \cdots ,|V|\}\f$ be a set of <EM>vertices</EM>. 
  A <EM>simplex</EM> \f$\sigma\f$ is a subset of vertices 
  \f$\sigma \subseteq V\f$. A <EM>simplicial complex</EM> \f$\mathbf{K}\f$ 
  on \f$V\f$ is a collection of simplices \f$\{\sigma\}\f$,
  \f$\sigma \subseteq V\f$, such that \f$\tau \subseteq \sigma \in \mathbf{K} 
  \Rightarrow \tau \in \mathbf{K}\f$. The dimension \f$n=|\sigma|-1\f$ of \f$\sigma\f$ 
  is its number of elements minus 1. A <EM>filtration</EM> of a simplicial complex is
  a function \f$f:\mathbf{K} \rightarrow \mathbb{R}\f$ satisfying \f$f(\tau)\leq 
  f(\sigma)\f$ whenever \f$\tau \subseteq \sigma\f$. 

  We define the concept FilteredComplex which enumerates the requirements for a class 
  to represent a filtered complex from which persistent homology may be computed.
    We use the vocabulary of simplicial complexes, but the concept 
    is valid for any type of cell complex. The main requirements 
    are the definition of:
    \li type <CODE>Indexing_tag</CODE>, which is a model of the concept 
        <CODE>IndexingTag</CODE>, 
        describing the nature of the indexing scheme,
    \li type Simplex_handle to manipulate simplices,
    \li method <CODE>int dimension(Simplex_handle)</CODE> returning 
        the dimension of a simplex,
    \li type and method <CODE>Boundary_simplex_range 
        boundary_simplex_range(Simplex_handle)</CODE> that returns 
        a range giving access to the codimension 1 subsimplices of the 
        input simplex, as-well-as the coefficients \f$(-1)^i\f$ in the 
        definition of the operator \f$\partial\f$. The iterators have 
        value type <CODE>Simplex_handle</CODE>,
    \li type and method 
        <CODE>Filtration_simplex_range filtration_simplex_range ()</CODE> 
        that returns a range giving 
        access to all the simplices of the complex read in the order 
        assigned by the indexing scheme,
    \li type and method 
        <CODE>Filtration_value filtration (Simplex_handle)</CODE> that returns the value of 
        the filtration on the simplex represented by the handle.

    <DT>Homology:</DT> 
    For a ring \f$\mathcal{R}\f$, the group of <EM>n-chains</EM>, 
    denoted \f$\mathbf{C}_n(\mathbf{K},\mathcal{R})\f$, of \f$\mathbf{K}\f$ is the 
    group of formal sums of
    n-simplices with \f$\mathcal{R}\f$ coefficients. The <EM>boundary operator</EM> is a 
    linear operator 
    \f$\partial_n: \mathbf{C}_n(\mathbf{K},\mathcal{R}) \rightarrow \mathbf{C}_{n-1}(\mathbf{K},\mathcal{R})\f$ 
    such that \f$\partial_n \sigma = \partial_n [v_0, \cdots , v_n] = 
    \sum_{i=0}^n (-1)^{i}[v_0,\cdots ,\widehat{v_i}, \cdots,v_n]\f$,
    where \f$\widehat{v_i}\f$ means \f$v_i\f$ is omitted from the list. The chain 
    groups form a sequence:

    \f[\cdots \ \ \mathbf{C}_n(\mathbf{K},\mathcal{R}) \xrightarrow{\ \partial_n\ } \mathbf{C}_{n-1}(\mathbf{K},\mathcal{R}) 
    \xrightarrow{\partial_{n-1}} \cdots \xrightarrow{\ \partial_2 \ }  
    \mathbf{C}_1(\mathbf{K},\mathcal{R}) \xrightarrow{\ \partial_1 \ }  \mathbf{C}_0(\mathbf{K},\mathcal{R}) \f]
    
    of finitely many groups \f$\mathbf{C}_n(\mathbf{K},\mathcal{R})\f$ and homomorphisms 
    \f$\partial_n\f$, indexed by the dimension \f$n \geq 0\f$.
    The boundary operators satisfy the property \f$\partial_n \circ \partial_{n+1}=0\f$ 
    for every \f$n > 0\f$ 
    and we define the homology groups: 
    
    \f[\mathbf{H}_n(\mathbf{K},\mathcal{R}) = \ker \partial_n / \mathrm{im} \  \partial_{n+1}\f]
    
    We refer to \cite Munkres-elementsalgtop1984 for an introduction to homology 
    theory and to \cite DBLP:books/daglib/0025666 for an introduction to persistent homology.
    
    <DT>Indexing Scheme:</DT>
    "Changing" a simplicial complex consists in applying a simplicial map.
    An <EM>indexing scheme</EM> is a directed graph together with a traversal 
    order, such that two 
    consecutive nodes in the graph are connected by an arrow (either forward or backward). 
    The nodes represent simplicial complexes and the directed edges simplicial maps.

    From the computational point of view, there are two types of indexing schemes of 
    interest 
    in persistent homology: <EM>linear</EM> ones 
    \f$\bullet \longrightarrow \bullet \longrightarrow \cdots \longrightarrow \bullet 
    \longrightarrow \bullet\f$
    in persistent homology \cite DBLP:journals/dcg/ZomorodianC05 ,
    and <EM>zigzag</EM> ones 
    \f$\bullet \longrightarrow \bullet \longleftarrow \cdots 
    \longrightarrow \bullet 
    \longleftarrow \bullet \f$ in zigzag persistent 
    homology \cite DBLP:journals/focm/CarlssonS10. 
    These indexing schemes have a natural left-to-right traversal order, and we 
    describe them with ranges and iterators.
    In the current release of the Gudhi library, only the linear case is implemented.

    In the following, we consider the case where the indexing scheme is induced 
    by a filtration. 
    Ordering the simplices 
    by increasing filtration values (breaking ties so as a simplex appears after 
    its subsimplices of same filtration value) provides an indexing scheme.



    <DT>Implementations:</DT>
    We use the <EM>Compressed Annotation Matrix</EM> of \cite DBLP:conf/esa/BoissonnatDM13 to 
    implement the 
    persistent cohomology algorithm of \cite DBLP:journals/dcg/SilvaMV11 
    and \cite DBLP:conf/compgeom/DeyFW14 for persistence in the class Persistent_cohomology. 

    The coefficient fields available as models of CoefficientField are Field_Zp
    for \f$\mathbb{Z}_p\f$ (for any prime p) and Multi_field for the multi-field persistence algorithm 
    -- computing persistence simultaneously in various coefficient fields -- described 
    in \cite boissonnat:hal-00922572.


\section Examples
    We provide several example files: run these examples with -h for details on their use, and read the README file.

\li <CODE>rips_persistence.cpp</CODE> computes the Rips complex of a point cloud and its persistence diagram.

\li <CODE>rips_multifield_persistence.cpp</CODE> computes the Rips complex of a point cloud and its persistence diagram 
with a family of field coefficients.

\li <CODE>performance_rips_persistence.cpp</CODE> provides timings for the construction of the Rips complex on a set of 
points sampling a Klein bottle in \f$\mathbb{R}^5\f$ with a simplex tree, its conversion to a 
Hasse diagram and the computation of persistent homology and multi-field persistent homology for the 
different representations.



 \author    Clément Maria
 \version   1.0
 \date      2014
 \copyright GNU General Public License v3.
 @{  
 */

/** \brief Computes the persistent cohomology of a filtered complex.
*
* The computation is implemented with a Compressed Annotation Matrix 
* (CAM)\cite DBLP:conf/esa/BoissonnatDM13, 
* and is adapted to the computation of Multi-Field Persistent Homology (MF) 
* \cite boissonnat:hal-00922572 .
*
* \implements PersistentHomology
*
*/ 
//Memory allocation policy: classic, use a mempool, etc.*/
 template < class FilteredComplex
          , class CoefficientField
          > //to do mem allocation policy: classic, mempool, etc.
 class Persistent_cohomology {
 public:

  typedef FilteredComplex                           Complex_ds;
  // Data attached to each simplex to interface with a Property Map.
  typedef typename Complex_ds::Simplex_key          Simplex_key;   
  typedef typename Complex_ds::Simplex_handle       Simplex_handle;
  typedef typename Complex_ds::Filtration_value     Filtration_value;
  typedef typename CoefficientField::Element        Arith_element;
// Compressed Annotation Matrix types:
  // Column type
  typedef Persistent_cohomology_column < Simplex_key
                                       , Arith_element 
                                       >            Column; // contains 1 set_hook
  // Cell type
  typedef typename Column::Cell                     Cell;   // contains 2 list_hooks
  // Remark: constant_time_size must be false because base_hook_cam_h has auto_unlink link_mode
  typedef boost::intrusive::list < Cell
                                 , boost::intrusive::constant_time_size<false> 
                                 , boost::intrusive::base_hook< base_hook_cam_h >    
                                 >                  Hcell;  

  typedef boost::intrusive::set < Column
                                , boost::intrusive::constant_time_size<false> 
                                >                   Cam;
// Sparse column type for the annotation of the boundary of an element.
  typedef std::vector< std::pair<Simplex_key
                                , Arith_element > > A_ds_type;
// Persistent interval type. The Arith_element field is used for the multi-field framework.
  typedef boost::tuple< Simplex_handle
                      , Simplex_handle
                      , Arith_element >             Persistent_interval;
 
/** \brief Initializes the Persistent_cohomology class.
  *
  * @param[in] cpx Complex for which the persistent homology is compiuted. 
                   cpx is a model of FilteredComplex
  *
  * @param[in] persistence_dim_max if true, the persistent homology for the maximal dimension in the 
  *                                complex is computed. If false, it is ignored. Default is false.
  */
Persistent_cohomology ( Complex_ds & cpx 
                      , bool         persistence_dim_max = false )
: cpx_    (&cpx)
, dim_max_(cpx.dimension())                   // upper bound on the dimension of the simplices
, coeff_field_()                              // initialize the field coefficient structure.
, ds_rank_  (cpx_->num_simplices())           // union-find
, ds_parent_(cpx_->num_simplices())           // union-find
, ds_repr_  (cpx_->num_simplices(),NULL)      // union-find -> annotation vectors
, dsets_(&ds_rank_[0],&ds_parent_[0])         // union-find
, cam_()                                      // collection of annotation vectors
, zero_cocycles_()                            // union-find -> Simplex_key of creator for 0-homology
, transverse_idx_()                           // key -> row
, persistent_pairs_() 
, interval_length_policy(&cpx,0)
, column_pool_(new boost::object_pool< Column > ()) // memory pools for the CAM
, cell_pool_(new boost::object_pool< Cell > ())
{
  if( persistence_dim_max ) { ++dim_max_; }
  Simplex_key idx_fil = 0;
  for(auto & sh : cpx_->filtration_simplex_range())
  { 
    cpx_->assign_key(sh,idx_fil); ++idx_fil;
    dsets_.make_set(cpx_->key(sh)); 
  }
}

~Persistent_cohomology()
{ 
//Clean the remaining columns in the matrix.
  for(auto & cam_ref : cam_) { cam_ref.col_.clear(); }
//Clean the transversal lists
  for( auto & transverse_ref : transverse_idx_ )
  {  transverse_ref.second.row_->clear(); delete transverse_ref.second.row_; }
//Clear the memory pools
  delete column_pool_;
  delete cell_pool_;
}

private:
struct length_interval {
  length_interval ( Complex_ds * cpx
                  , Filtration_value min_length)
  : cpx_(cpx)
  , min_length_(min_length) {}

  bool operator()(Simplex_handle sh1, Simplex_handle sh2)
  { return cpx_->filtration(sh2) - cpx_->filtration(sh1) > min_length_; }

  void set_length(Filtration_value new_length) { min_length_ = new_length; }

  Complex_ds          * cpx_;
  Filtration_value      min_length_;
};


public:
/** \brief Initializes the coefficient field.*/
void init_coefficients( int charac     ) { coeff_field_.init(charac); }
/** \brief Initializes the coefficient field for multi-field persistent homology.*/
void init_coefficients( int charac_min
                      , int charac_max ) { coeff_field_.init(charac_min,charac_max); }

/** \brief Compute the persistent homology of the filtered simplicial
  * complex.
  *
  * @param[in] min_interval_length the computation disgards all intervals of length
  *                                less or equal than min_interval_length
  *
  * Assumes that the filtration provided by the simplicial complex is 
  * valid. Undefined behavior otherwise. */
void compute_persistent_cohomology ( Filtration_value min_interval_length = 0 )
{
  interval_length_policy.set_length(min_interval_length);
  // Compute all finite intervals
  for( auto sh : cpx_->filtration_simplex_range() )
  {
    int dim_simplex = cpx_->dimension(sh);
    switch(dim_simplex) {
      case 0 :                                              break;
      case 1 : update_cohomology_groups_edge( sh )        ; break;
      default: update_cohomology_groups( sh, dim_simplex ); break;
    }
  }
  // Compute infinite intervals of dimension 0
  Simplex_key key;
  for(auto v_sh : cpx_->skeleton_simplex_range(0)) //for all 0-dimensional simplices
  { 
    key = cpx_->key(v_sh);

    if( ds_parent_[key] == key  //root of its tree
        && zero_cocycles_.find(key) == zero_cocycles_.end() ) 
    {
      persistent_pairs_.push_back( Persistent_interval ( cpx_->simplex(key)
                                                       , cpx_->null_simplex()
                                                       , coeff_field_.characteristic() )
                                 );
    }
  }
  for( auto zero_idx : zero_cocycles_ )
  {
    persistent_pairs_.push_back( Persistent_interval ( cpx_->simplex(zero_idx.second)
                                                     , cpx_->null_simplex()
                                                     , coeff_field_.characteristic() )
                               );
  }
// Compute infinite interval of dimension > 0  
  for(auto cocycle : transverse_idx_)
  {
    persistent_pairs_.push_back( Persistent_interval ( cpx_->simplex (cocycle.first)
                                                     , cpx_->null_simplex()
                                                     , cocycle.second.characteristics_ ) );
  }
}



private:
/** \brief Update the cohomology groups under the insertion of an edge.
  * 
  * The 0-homology is maintained with a simple Union-Find data structure, which
  * explains the existance of a specific function of edge insertions. */
void update_cohomology_groups_edge ( Simplex_handle sigma ) 
{
  Simplex_handle u,v;
  boost::tie(u,v) = cpx_->endpoints(sigma);
  
  Simplex_key ku = dsets_.find_set( cpx_->key(u) ); 
  Simplex_key kv = dsets_.find_set( cpx_->key(v) );

  if(ku != kv ) {        // Destroy a connected component
    dsets_.link(ku,kv);        
    // Keys of the simplices which created the connected components containing
    // respectively u and v. 
    Simplex_key idx_coc_u, idx_coc_v;
    auto map_it_u = zero_cocycles_.find(ku);
    // If the index of the cocycle representing the class is already ku.
    if (map_it_u == zero_cocycles_.end()) { idx_coc_u = ku;               }
    else                                  { idx_coc_u = map_it_u->second; }

    auto map_it_v = zero_cocycles_.find(kv);
    // If the index of the cocycle representing the class is already kv.
    if (map_it_v == zero_cocycles_.end()) { idx_coc_v = kv;               }
    else                                  { idx_coc_v = map_it_v->second; }

    if(cpx_->filtration(cpx_->simplex(idx_coc_u)) 
        < cpx_->filtration(cpx_->simplex(idx_coc_v)) ) // Kill cocycle [idx_coc_v], which is younger. 
      {
        if(interval_length_policy(cpx_->simplex(idx_coc_v),sigma)) {
          persistent_pairs_.push_back ( Persistent_interval ( cpx_->simplex(idx_coc_v)
                                                            , sigma
                                                            , coeff_field_.characteristic() 
                                                            )
                                      );
      }
    // Maintain the index of the 0-cocycle alive.
        if( kv != idx_coc_v ) { zero_cocycles_.erase( map_it_v ); }
        if( kv == dsets_.find_set(kv) ) {
          if( ku != idx_coc_u ) { zero_cocycles_.erase( map_it_u ); }
          zero_cocycles_[kv] = idx_coc_u;
        }
      }
    else // Kill cocycle [idx_coc_u], which is younger.
      {
        if(interval_length_policy(cpx_->simplex(idx_coc_u),sigma)) {
          persistent_pairs_.push_back ( Persistent_interval ( cpx_->simplex(idx_coc_u)
                                                            , sigma
                                                            , coeff_field_.characteristic() 
                                                            )
                                      );
        }
    // Maintain the index of the 0-cocycle alive.
        if( ku != idx_coc_u ) { zero_cocycles_.erase( map_it_u ); }
        if( ku == dsets_.find_set(ku) ) {
          if( kv != idx_coc_v ) { zero_cocycles_.erase( map_it_v ); }
          zero_cocycles_[ku] = idx_coc_v;
        }
      }
    cpx_->assign_key(sigma,cpx_->null_key()); 
  }
  else { // If ku == kv, same connected component: create a 1-cocycle class.
    create_cocycle( sigma, coeff_field_.multiplicative_identity(), coeff_field_.characteristic() ); 
  } 
}

/*
 * Compute the annotation of the boundary of a simplex.
 */
void annotation_of_the_boundary(std::map< Simplex_key, Arith_element > & map_a_ds
                               , Simplex_handle sigma
                               , int dim_sigma )
{
  // traverses the boundary of sigma, keeps track of the annotation vectors,
  // with multiplicity, in a map.
  std::map < Column *, int >                  annotations_in_boundary;
  std::pair < typename std::map< Column *, int >::iterator
            , bool >                          result_insert_bound;
  int sign = 1 - 2 * (dim_sigma % 2); // \in {-1,1} provides the sign in the 
                                      // alternate sum in the boundary.
  Simplex_key key;      Column * curr_col;

  for( auto sh : cpx_->boundary_simplex_range(sigma) )
  {
    key = cpx_->key(sh);
    if( key != cpx_->null_key() ) // A simplex with null_key is a killer, and have null annotation
      {                           // vector.
        // Find its annotation vector
        curr_col = ds_repr_[ dsets_.find_set(key) ];
        if( curr_col != NULL ) 
        { // and insert it in annotations_in_boundary with multyiplicative factor "sign".
          result_insert_bound = 
            annotations_in_boundary.insert(std::pair<Column *,int>(curr_col,sign));  
          if( !(result_insert_bound.second) ) { result_insert_bound.first->second += sign; }
        }
      }
    sign = -sign;
  }
  // Sum the annotations with multiplicity, using a map<key,coeff> 
  // to represent a sparse vector.
  std::pair < typename std::map < Simplex_key, Arith_element >::iterator
            , bool >                                                result_insert_a_ds;

  for( auto ann_ref : annotations_in_boundary ) 
  {    
    if(ann_ref.second != coeff_field_.additive_identity()) // For all columns in the boundary,
    {                                            
      for( auto cell_ref : ann_ref.first->col_ ) // insert every cell in map_a_ds with multiplicity
      { 
        Arith_element w_y = 
            coeff_field_.times(cell_ref.coefficient_ , ann_ref.second); //coefficient * multiplicity

        if( w_y != coeff_field_.additive_identity() ) // if != 0
        {
          result_insert_a_ds = map_a_ds.insert(std::pair< Simplex_key
                                                        , Arith_element >(cell_ref.key_ , w_y));
          if( !(result_insert_a_ds.second) )   //if cell_ref.key_ already a Key in map_a_ds
          { 
            coeff_field_.plus_equal(result_insert_a_ds.first->second, w_y); 
            if(result_insert_a_ds.first->second == coeff_field_.additive_identity())
              {  map_a_ds.erase(result_insert_a_ds.first); }
          }
        }
      }
    }  
  }
}

/* 
 * Update the cohomology groups under the insertion of a simplex.
 */
void update_cohomology_groups ( Simplex_handle sigma
                              , int dim_sigma )
{
//Compute the annotation of the boundary of sigma:
  std::map< Simplex_key, Arith_element >            map_a_ds;
  annotation_of_the_boundary(map_a_ds, sigma, dim_sigma );
// Update the cohomology groups:
  if( map_a_ds.empty() ) {  // sigma is a creator in all fields represented in coeff_field_
    if(dim_sigma < dim_max_) { create_cocycle ( sigma
                                              , coeff_field_.multiplicative_identity() 
                                              , coeff_field_.characteristic() );}
  }
  else {                    // sigma is a destructor in at least a field in coeff_field_
 // Convert map_a_ds to a vector
    A_ds_type a_ds; //admits reverse iterators
    for ( auto map_a_ds_ref : map_a_ds )
    { 
      a_ds.push_back( std::pair< Simplex_key
                               , Arith_element> ( map_a_ds_ref.first
                                                , map_a_ds_ref.second ));
    }

   
    Arith_element inv_x, charac;
    Arith_element prod = coeff_field_.characteristic(); // Product of characteristic of the fields
    for( auto a_ds_rit = a_ds.rbegin(); 
         (a_ds_rit != a_ds.rend()) && (prod != coeff_field_.multiplicative_identity());
         ++a_ds_rit )
    {
      std::tie(inv_x,charac) = coeff_field_.inverse ( a_ds_rit->second
                                                    , prod );   

      if( inv_x != coeff_field_.additive_identity() )
      {
        destroy_cocycle ( sigma
                        , a_ds
                        , a_ds_rit->first
                        , inv_x 
                        , charac );
        prod /= charac;
      }
    }
    if( prod != coeff_field_.multiplicative_identity() && dim_sigma < dim_max_ )
      { create_cocycle( sigma , coeff_field_.multiplicative_identity(prod), prod ); }
  }
}

/*  \brief Create a new cocycle class.
  *
  * The class is created by the insertion of the simplex sigma.
  * The methods adds a cocycle, representing the new cocycle class,
  * to the matrix representing the cohomology groups.
  * The new cocycle has value 0 on every simplex except on sigma
  * where it worths 1.*/
void create_cocycle ( Simplex_handle sigma
                    , Arith_element  x 
                    , Arith_element  charac )
{
  Simplex_key key = cpx_->key(sigma);
  // Create a column containing only one cell,
  Column * new_col  = column_pool_->construct(Column(key)); 
  Cell   * new_cell = cell_pool_->construct(Cell (key, x, new_col)); 
  new_col->col_.push_back(*new_cell);  
  // and insert it in the matrix, in constant time thanks to the hint cam_.end().
  // Indeed *new_col has the biggest lexicographic value because key is the 
  // biggest key used so far.
  cam_.insert (cam_.end(), *new_col); 
  // Update the disjoint sets data structure.
  Hcell * new_hcell = new Hcell;
  new_hcell->push_back(*new_cell);
  transverse_idx_[key] = cocycle(charac,new_hcell); //insert the new row
  ds_repr_[key] = new_col;
}

/*  \brief Destroy a cocycle class.
  *
  * The cocycle class is destroyed by the insertion of sigma.
  * The methods proceeds to a reduction of the matrix representing 
  * the cohomology groups using Gauss pivoting. The reduction zeros-out
  * the row containing the cell with highest key in
  * a_ds, the annotation of the boundary of simplex sigma. This key
  * is "death_key".*/
void destroy_cocycle ( Simplex_handle   sigma
                     , A_ds_type const& a_ds 
                     , Simplex_key      death_key
                     , Arith_element    inv_x
                     , Arith_element    charac )
{
  // Create a finite persistent interval
  if(interval_length_policy(cpx_->simplex(death_key),sigma)) {
    persistent_pairs_.push_back ( Persistent_interval ( cpx_->simplex(death_key) //creator
                                                      , sigma                    //destructor
                                                      , charac )                 //fields 
                                );                              // for which the interval exists 
  }

  auto death_key_row = transverse_idx_.find(death_key); // Find the beginning of the row.
  std::pair< typename Cam::iterator, bool > result_insert_cam;

  auto row_cell_it = death_key_row->second.row_->begin();
  
  while( row_cell_it != death_key_row->second.row_->end() ) // Traverse all cells in 
  {                                                         // the row at index death_key.
    Arith_element w = coeff_field_.times_minus( inv_x , row_cell_it->coefficient_ );

    if( w != coeff_field_.additive_identity() ) 
    { 
      Column * curr_col = row_cell_it->self_col_;         ++row_cell_it;
      // Disconnect the column from the rows in the CAM.
      for( auto col_cell_it = curr_col->col_.begin();
          col_cell_it != curr_col->col_.end(); ++col_cell_it ) 
        { col_cell_it->base_hook_cam_h::unlink(); }

      // Remove the column from the CAM before modifying its value
      cam_.erase( cam_.iterator_to(*curr_col) ); 
      // Proceed to the reduction of the column
      plus_equal_column(*curr_col, a_ds, w);

      if( curr_col->col_.empty() ) // If the column is null
      { 
        ds_repr_[ curr_col->class_key_ ] = NULL;  
        column_pool_->free(curr_col); //delete curr_col; 
      }
      else 
      { 
        // Find whether the column obtained is already in the CAM
        result_insert_cam = cam_.insert( *curr_col );
        if ( result_insert_cam.second ) // If it was not in the CAM before: insertion has succeeded
        {
          for ( auto col_cell_it = curr_col->col_.begin();
                col_cell_it != curr_col->col_.end(); ++col_cell_it ) //re-establish the row links
          { transverse_idx_[ col_cell_it->key_ ].row_->push_front(*col_cell_it); }
        }
        else // There is already an identical column in the CAM: 
        {    // merge two disjoint sets.
          dsets_.link ( curr_col->class_key_ , 
                        result_insert_cam.first->class_key_ );

          Simplex_key key_tmp = dsets_.find_set( curr_col->class_key_ );
          ds_repr_[ key_tmp ] = &(*(result_insert_cam.first));
          result_insert_cam.first->class_key_ = key_tmp;
          column_pool_->free(curr_col); //delete curr_col;
        }
      }
    }
    else { ++row_cell_it; } // If w == 0, pass.
  }

  // Because it is a killer simplex, set the data of sigma to null_key().
  if(charac == coeff_field_.characteristic()) { cpx_->assign_key( sigma, cpx_->null_key() ); }
  if(death_key_row->second.characteristics_ == charac) 
  { 
    delete death_key_row->second.row_;
    transverse_idx_.erase(death_key_row); 
  }
  else { death_key_row->second.characteristics_ /= charac; }
}

/* 
 * Assign:    target <- target + w * other. 
 */
void plus_equal_column ( Column & target
                       , A_ds_type const& other //value_type is pair<Simplex_key,Arith_element>
                       , Arith_element w )
{
  auto target_it = target.col_.begin(); auto other_it = other.begin();
  while ( target_it != target.col_.end() && other_it != other.end() )
  {
    if(target_it->key_ < other_it->first) { ++target_it; }
    else {
      if(target_it->key_ > other_it->first) 
      {
        Cell * cell_tmp = cell_pool_->construct(Cell( other_it->first   //key
                                                    , coeff_field_.additive_identity()
                                                    , &target));

        coeff_field_.plus_times_equal(cell_tmp->coefficient_, other_it->second, w);

        target.col_.insert( target_it, *cell_tmp );

        ++other_it;
      }
      else { //it1->key == it2->key
        //target_it->coefficient_ <- target_it->coefficient_ + other_it->second * w
        coeff_field_.plus_times_equal( target_it->coefficient_ , other_it->second , w);
        if( target_it->coefficient_ == coeff_field_.additive_identity() )
        {
          auto tmp_it = target_it;
          ++target_it; ++other_it;   // iterators remain valid
          Cell * tmp_cell_ptr = &(*tmp_it); 
          target.col_.erase(tmp_it); // removed from column
        
          coeff_field_.clear_coefficient(tmp_cell_ptr->coefficient_);
          cell_pool_->free(tmp_cell_ptr); // delete from memory 
        }
        else { ++target_it; ++other_it; }
      }
    }
  }
  while(other_it != other.end()) 
  {
    Cell * cell_tmp = cell_pool_->construct(Cell( other_it->first   //key
                                                , coeff_field_.additive_identity()
                                                , &target));

    coeff_field_.plus_times_equal(cell_tmp->coefficient_, other_it->second, w);

    target.col_.insert( target.col_.end(), *cell_tmp );

    ++other_it;
  }
}

/*
 * Compare two intervals by length.
 */
struct cmp_intervals_by_length {
  cmp_intervals_by_length( Complex_ds * sc ) : sc_ (sc) {}
  bool operator() (  Persistent_interval & p1
                  ,  Persistent_interval & p2 )
  {
    return ( sc_->filtration( get<1>(p1) ) - sc_->filtration( get<0>(p1) ) 
             > sc_->filtration( get<1>(p2) ) - sc_->filtration( get<0>(p2) ) );
  }
  Complex_ds * sc_;
};

public:
/** \brief Output the persistence diagram in ostream.
  *
  * The file format is the following:
  *    p1*...*pr   dim b d 
  *
  * where "dim" is the dimension of the homological feature,
  * b and d are respectively the birth and death of the feature and
  * p1*...*pr is the product of prime numbers pi such that the homology 
  * feature exists in homology with Z/piZ coefficients.
  */
void output_diagram(std::ostream& ostream = std::cout)
{
  cmp_intervals_by_length cmp( cpx_ );
  persistent_pairs_.sort( cmp );
  for(auto pair : persistent_pairs_)
  {
    ostream << get<2>(pair)                   << "  "
            << cpx_->dimension(get<0>(pair))  << " "
            << cpx_->filtration(get<0>(pair)) << " " 
            << cpx_->filtration(get<1>(pair)) << " " 
            << std::endl;
  }
}

private:
/* 
 * Structure representing a cocycle.
 */
struct cocycle {
  cocycle() {}
  cocycle( Arith_element characteristics
         , Hcell       * row )
  : row_(row), characteristics_(characteristics) {}

  Hcell * row_;                   //points to the corresponding row in the CAM
  Arith_element characteristics_; //product of field characteristics for which the cocycle exist
};

public:
  Complex_ds *         cpx_;
  int                  dim_max_;
  CoefficientField     coeff_field_;

/*  Disjoint sets data structure to link the model of FilteredComplex
  * with the compressed annotation matrix.
  * ds_rank_ is a property map Simplex_key -> int, ds_parent_ is a property map 
  * Simplex_key -> simplex_key_t */  
  std::vector< int >                            ds_rank_;  
  std::vector< Simplex_key >                    ds_parent_;
  std::vector< Column * >                       ds_repr_;
  boost::disjoint_sets< int *, Simplex_key * >  dsets_;
/* The compressed annotation matrix fields.*/
  Cam                                           cam_;
/*  Dictionary establishing the correspondance between the Simplex_key of
  * the root vertex in the union-find ds and the Simplex_key of the vertex which
  * created the connected component as a 0-dimension homology feature.*/
  std::map<Simplex_key,Simplex_key>             zero_cocycles_;
/*  Key -> row. */ 
  std::map< Simplex_key , cocycle >             transverse_idx_;
/* Persistent intervals. */
  std::list< Persistent_interval >              persistent_pairs_; 
  length_interval                               interval_length_policy;

  boost::object_pool< Column > *                column_pool_;
  boost::object_pool< Cell >   *                cell_pool_;
};

/** @} */ //end defgroup persistent_cohomology

}  // namespace GUDHI

#endif // _PERSISTENCECOMPUTATION_SIMPLEXTREE_
