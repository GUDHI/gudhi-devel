/*
 * Skeleton_blocker_link_superior.h
 *
 *  Created on: Feb 19, 2014 
 *      Author: David Salinas
 *  Copyright 2013 INRIA. All rights reserved
 */

#ifndef GUDHI_SKELETON_BLOCKER_LINK_SUPERIOR_H_
#define GUDHI_SKELETON_BLOCKER_LINK_SUPERIOR_H_

#include "gudhi/Skeleton_blocker_link_complex.h"

namespace Gudhi{

namespace skbl {

template<class ComplexType> class Skeleton_blocker_sub_complex;


/**
 *  \brief Class representing the link of a simplicial complex encoded by a skeleton/blockers pair.
 *  It computes only vertices greater than the simplex used to build the link.
 */
template<typename ComplexType>
class Skeleton_blocker_link_superior : public Skeleton_blocker_link_complex<ComplexType>
{
	typedef typename ComplexType::Edge_handle Edge_handle;

	typedef typename ComplexType::boost_vertex_handle boost_vertex_handle;

public:

	typedef typename ComplexType::Vertex_handle Vertex_handle;
	typedef typename ComplexType::Root_vertex_handle Root_vertex_handle;
	typedef typename ComplexType::Simplex_handle Simplex_handle;
	typedef typename ComplexType::Root_simplex_handle Root_simplex_handle;
	typedef typename ComplexType::BlockerMap BlockerMap;
	typedef typename ComplexType::BlockerPair BlockerPair;
	typedef typename ComplexType::BlockerMapIterator BlockerMapIterator;
	typedef typename ComplexType::BlockerMapConstIterator BlockerMapConstIterator;
	typedef typename ComplexType::Simplex_handle::Simplex_vertex_const_iterator AddressSimplexConstIterator;
	typedef typename ComplexType::Root_simplex_handle::Simplex_vertex_const_iterator IdSimplexConstIterator;


	Skeleton_blocker_link_superior()
		:Skeleton_blocker_link_complex<ComplexType>(true)
		 {
		 }

	Skeleton_blocker_link_superior(const ComplexType & parent_complex, Simplex_handle& alpha_parent_adress)
		:Skeleton_blocker_link_complex<ComplexType>(parent_complex,alpha_parent_adress,true)
		 {
		 }

	Skeleton_blocker_link_superior(const ComplexType & parent_complex, Vertex_handle a_parent_adress)
	:Skeleton_blocker_link_complex<ComplexType>(parent_complex,a_parent_adress,true)
	 {
	}

};

}

}  // namespace GUDHI


#endif /* SKELETON_BLOCKER_LINK_SUPERIOR_H_ */
