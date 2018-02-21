/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2015  INRIA Sophia Antipolis-Méditerranée (France)
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

#ifndef HASSE_DIAGRAM_H_
#define HASSE_DIAGRAM_H_

#include <vector>
#include <list>

class Hasse_diagram {

private:

  struct Node {
    typedef typename std::list<Node>::iterator Node_iterator; 
    std::list<Node_iterator> facets_, cofaces_;
  };
    
  typedef std::vector<std::list<Node> > Dim_lists;

  Dim_lists dim_lists_;

  // debug
  int max_simplices_out = 0;
  
public:

  typedef typename std::list<Node>::iterator Node_iterator;

  Hasse_diagram() {    
  }

  Node_iterator add_node() {
    if (dim_lists_.empty())
      dim_lists_.push_back(std::list<Node>());
    dim_lists_[0].emplace_back(Node());
    return dim_lists_[0].end()--;
  }
  
}; //class Hasse_diagram


#endif
