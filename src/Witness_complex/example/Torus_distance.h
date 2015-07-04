#ifndef GUDHI_TORUS_DISTANCE_H_
#define GUDHI_TORUS_DISTANCE_H_

#include <math.h>

#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Epick_d.h>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::Point_d Point_d;
typedef K::FT FT;
typedef CGAL::Search_traits<
  FT, Point_d,
  typename K::Cartesian_const_iterator_d,
  typename K::Construct_cartesian_const_iterator_d> Traits_base;

/**
 * \brief Class of distance in a flat torus in dimension D
 *
 */
class Torus_distance {

public:
    typedef K::FT    FT;
    typedef K::Point_d Point_d;
    typedef Point_d Query_item;
    typedef typename CGAL::Dynamic_dimension_tag D;

  double box_length = 2;

  FT transformed_distance(Query_item q, Point_d p) const
  {
    FT distance = FT(0);
    FT coord = FT(0);
    //std::cout << "Hello skitty!\n";
    typename K::Construct_cartesian_const_iterator_d construct_it=Traits_base().construct_cartesian_const_iterator_d_object();
    typename K::Cartesian_const_iterator_d qit = construct_it(q),
	qe = construct_it(q,1), pit = construct_it(p);
	for(; qit != qe; qit++, pit++)
          {
            coord = sqrt(((*qit)-(*pit))*((*qit)-(*pit)));
            if (coord*coord <= (box_length-coord)*(box_length-coord))
              distance += coord*coord;
            else
              distance += (box_length-coord)*(box_length-coord);
          }
        return distance;
  }

  FT min_distance_to_rectangle(const Query_item& q,
                               const CGAL::Kd_tree_rectangle<FT,D>& r) const {
    FT distance = FT(0);
    FT dist1, dist2;
    typename K::Construct_cartesian_const_iterator_d construct_it=Traits_base().construct_cartesian_const_iterator_d_object();
    typename K::Cartesian_const_iterator_d qit = construct_it(q),
                                           qe = construct_it(q,1);
    for(unsigned int i = 0;qit != qe; i++, qit++)
      {
        if((*qit) < r.min_coord(i))
          {
            dist1 = (r.min_coord(i)-(*qit));
            dist2 = (box_length - r.max_coord(i)+(*qit));
            if (dist1 < dist2)
              distance += dist1*dist1;
            else
              distance += dist2*dist2;
          }
        else if ((*qit) > r.max_coord(i))
          {
            dist1 = (box_length - (*qit)+r.min_coord(i));
            dist2 = ((*qit) - r.max_coord(i));
            if (dist1 < dist2)
              distance += dist1*dist1;
            else
              distance += dist2*dist2;
          }
      }
    return distance;
  }

  FT min_distance_to_rectangle(const Query_item& q,
                               const CGAL::Kd_tree_rectangle<FT,D>& r,
                               std::vector<FT>& dists) const {
    FT distance = FT(0);
    FT dist1, dist2;
    typename K::Construct_cartesian_const_iterator_d construct_it=Traits_base().construct_cartesian_const_iterator_d_object();
    typename K::Cartesian_const_iterator_d qit = construct_it(q),
	                  	           qe = construct_it(q,1);
    //std::cout << r.max_coord(0) << std::endl;
    for(unsigned int i = 0;qit != qe; i++, qit++)
      {
        if((*qit) < r.min_coord(i))
          {
            dist1 = (r.min_coord(i)-(*qit));
            dist2 = (box_length - r.max_coord(i)+(*qit));
            if (dist1 < dist2)
              {
                dists[i] = dist1;
                distance += dist1*dist1;
              }
            else
              {
                dists[i] = dist2;
                distance += dist2*dist2;
                //std::cout << "Good stuff1\n";
              }
          }
        else if ((*qit) > r.max_coord(i))
          {
            dist1 = (box_length - (*qit)+r.min_coord(i));
            dist2 = ((*qit) - r.max_coord(i));
            if (dist1 < dist2)
              {
                dists[i] = dist1;
                distance += dist1*dist1;
                //std::cout << "Good stuff2\n";
              }
            else
              {
                dists[i] = dist2;
                distance += dist2*dist2;
              }
          }
      };
    return distance;
  }
    
  FT max_distance_to_rectangle(const Query_item& q,
                               const CGAL::Kd_tree_rectangle<FT,D>& r) const {
    FT distance=FT(0);
    typename K::Construct_cartesian_const_iterator_d construct_it=Traits_base().construct_cartesian_const_iterator_d_object();
    typename K::Cartesian_const_iterator_d qit = construct_it(q),
		                           qe = construct_it(q,1);
    for(unsigned int i = 0;qit != qe; i++, qit++)
      {
        if (box_length <= (r.min_coord(i)+r.max_coord(i)))
          if ((r.max_coord(i)+r.min_coord(i)-box_length)/FT(2.0) <= (*qit) &&
              (*qit) <= (r.min_coord(i)+r.max_coord(i))/FT(2.0))
            distance += (r.max_coord(i)-(*qit))*(r.max_coord(i)-(*qit));
          else
            distance += ((*qit)-r.min_coord(i))*((*qit)-r.min_coord(i));
        else
          if ((box_length-r.max_coord(i)-r.min_coord(i))/FT(2.0) <= (*qit) ||
              (*qit) <= (r.min_coord(i)+r.max_coord(i))/FT(2.0))
            distance += (r.max_coord(i)-(*qit))*(r.max_coord(i)-(*qit));
          else
            distance += ((*qit)-r.min_coord(i))*((*qit)-r.min_coord(i));
      }
    return distance;
  }

  
  FT max_distance_to_rectangle(const Query_item& q,
                               const CGAL::Kd_tree_rectangle<FT,D>& r,
                               std::vector<FT>& dists) const {
    FT distance=FT(0);
    typename K::Construct_cartesian_const_iterator_d construct_it=Traits_base().construct_cartesian_const_iterator_d_object();
    typename K::Cartesian_const_iterator_d qit = construct_it(q),
		                           qe = construct_it(q,1);
    for(unsigned int i = 0;qit != qe; i++, qit++)
      {
        if (box_length <= (r.min_coord(i)+r.max_coord(i)))
          if ((r.max_coord(i)+r.min_coord(i)-box_length)/FT(2.0) <= (*qit) &&
              (*qit) <= (r.min_coord(i)+r.max_coord(i))/FT(2.0))
            {
              dists[i] = r.max_coord(i)-(*qit);
              distance += (r.max_coord(i)-(*qit))*(r.max_coord(i)-(*qit));
            }
          else
            {
              dists[i] = sqrt(((*qit)-r.min_coord(i))*((*qit)-r.min_coord(i)));
              distance += ((*qit)-r.min_coord(i))*((*qit)-r.min_coord(i));
            }
        else
          if ((box_length-r.max_coord(i)-r.min_coord(i))/FT(2.0) <= (*qit) ||
              (*qit) <= (r.min_coord(i)+r.max_coord(i))/FT(2.0))
            {
              dists[i] = sqrt((r.max_coord(i)-(*qit))*(r.max_coord(i)-(*qit)));
              distance += (r.max_coord(i)-(*qit))*(r.max_coord(i)-(*qit));
              
            }
          else
            {
              dists[i] = (*qit)-r.min_coord(i);
              distance += ((*qit)-r.min_coord(i))*((*qit)-r.min_coord(i));
            }
      }
    return distance;
  }

    inline FT new_distance(FT dist, FT old_off, FT new_off,
                           int )  const {
      
      FT new_dist = dist + (new_off*new_off - old_off*old_off);
      return new_dist;
    }
    
  inline FT transformed_distance(FT d) const {
    return d*d;
  }

  inline FT inverse_of_transformed_distance(FT d) const {
    return sqrt(d);
  }
  
};

#endif
