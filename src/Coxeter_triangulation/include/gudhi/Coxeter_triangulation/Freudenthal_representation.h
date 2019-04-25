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

  bool operator==(const Freudenthal_representation& other) const {
    if (dimension() != other.dimension())
      return false;
    if (vertex != other.vertex)
      return false;
    for (std::size_t k = 0; k < partition.size(); ++k)
      if (partition[k] != other.partition[k])
	return false;
    return true;
  }
};



#endif
