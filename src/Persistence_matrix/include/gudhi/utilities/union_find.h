/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef UNIONFIND_H
#define UNIONFIND_H

/**
 * @file Union_find.h
 * Contains the declaration of the \ref Union_find class.
 **/

#include <vector>
#include <unordered_map>
#include <iostream>

namespace Gudhi {
namespace persistence_matrix {

/**
 * class Union_find
 * @brief A basic Union-Find structure.
 * @tparam Element Type of the elements stored in the structure.
 * Has to provide a corresponding std::hash function.
 */
class Union_find {
public:
	using index = unsigned int;
	/**
	 * @brief Constructor.
	 */
	Union_find();
	Union_find(unsigned int numberOfElements);
	Union_find(const Union_find& toCopy);
	Union_find(Union_find&& other) noexcept;

	/**
	 * @brief Adds an element to the structure in its own component,
	 * when not already present.
	 * @param e Element to add.
	 */
	void initialize(unsigned int id);

	/**
	 * @brief Merges two components together
	 * if the first given element is contained.
	 * @param e1 An element of the first component.
	 * @param e2 An element of the second component.
	 * @pre e1 is contained in the structure.
	 */
	int merge(unsigned int id1, unsigned int id2);
	int find_and_merge(unsigned int id1, unsigned int id2);

	/**
	 * @brief Returns the representative of the component
	 * the given element is contained in.
	 * @pre @p e has to be contained in the structure.
	 * @param e Element to find.
	 * @return The representative of the component @p e is contained in.
	 */
	unsigned int find(unsigned int id);

	friend void swap(Union_find& uf1, Union_find& uf2){
		uf1.nodes_.swap(uf2.nodes_);
	}

private:
	struct Node{
		Node()
			: parent_(0), rank_(1){}
		Node(index parent, int rank)
			: parent_(parent), rank_(rank){}
		Node(const Node& toCopy)
			: parent_(toCopy.parent_), rank_(toCopy.rank_){}

//		Node& operator=(const Node& other){
//			parent_ = other.parent_;
//			rank_ = other.rank_;
//			return *this;
//		}

		index parent_;      ///< Position of the parent.
		int rank_;          ///< Rank in the parent hierachy.
	};

	std::vector<Node> nodes_;                   ///< Node container.
};

Union_find::Union_find()
{}

Union_find::Union_find(unsigned int numberOfElements)
	: nodes_(numberOfElements)
{
	for (unsigned int i = 0; i < nodes_.size(); ++i){
		nodes_[i].parent_ = i;
	}
}

Union_find::Union_find(const Union_find& toCopy)
	: nodes_(toCopy.nodes_)
{}

Union_find::Union_find(Union_find&& other) noexcept
	: nodes_(std::move(other.nodes_))
{}

void Union_find::initialize(unsigned int id) {
	if (nodes_.size() <= id) {
		unsigned int start = nodes_.size();
		nodes_.resize(id + 1);
		for (unsigned int i = start; i < nodes_.size(); ++i){
			nodes_[i].parent_ = i;
		}
	}
}

int Union_find::merge(unsigned int id1, unsigned int id2) {
	Node& n1 = nodes_[id1];
	Node& n2 = nodes_[id2];

	if (n1.parent_ == n2.parent_) {
		return -1;
	}

	if (n1.rank_ < n2.rank_){
		n1.parent_ = n2.parent_;
	} else if (n1.rank_ > n2.rank_){
		n2.parent_ = n1.parent_;
	} else {
		n2.parent_ = n1.parent_;
		(n1.rank_)++;
	}

	return n1.parent_;
}

int Union_find::find_and_merge(unsigned int id1, unsigned int id2) {
	return merge(find(id1), find(id2));
}

unsigned int Union_find::find(unsigned int id) {
	Node &n = nodes_[id];
	if (n.parent_ != id)
		n.parent_ = find(n.parent_);
	return n.parent_;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // UNIONFIND_H
