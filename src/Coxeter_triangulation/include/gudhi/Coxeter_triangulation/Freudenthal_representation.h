#ifndef FREUDENTHAL_REPRESENTATION_H_
#define FREUDENTHAL_REPRESENTATION_H_

template <class Vertex,
	  class Ordered_set_partition>
struct Freudenthal_representation {
  typedef Vertex vertex_type;
  typedef Ordered_set_partition partition_type;
  
  Vertex vertex;
  Ordered_set_partition partition;

  unsigned dimension() const {
    return partition.size() - 1;
  }
};

#endif
