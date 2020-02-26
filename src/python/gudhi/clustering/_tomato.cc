// g++ -O3 -Wall -shared -fPIC `python3 -m pybind11 --includes` XXX.cpp -o tomato`python3-config --extension-suffix`
#include <boost/container/flat_map.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/transform_value_property_map.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <vector>
#include <unordered_map>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>

namespace py=pybind11;

template<class T, class=std::enable_if_t<std::is_integral<T>::value>>
int getint(int i){return i;}
template<class T, class=decltype(std::declval<T>().template cast<int>())>
int getint(T i){return i.template cast<int>();}

// Raw clusters are clusters obtained through single-linkage, no merging.

typedef int Point_index;
typedef int Cluster_index;
struct Merge { Cluster_index first, second; double persist; };

template<class Neighbors, class Density, class Order, class ROrder>
auto tomato(Point_index num_points, Neighbors const&neighbors, Density const& density, Order const& order, ROrder const& rorder){
  // point index --> index of raw cluster it belongs to
  std::vector<Cluster_index> raw_cluster; raw_cluster.reserve(num_points);
  // cluster index --> index of top point in the cluster
  Cluster_index n_raw_clusters = 0; // current number of raw clusters seen
  //
  std::vector<Merge> merges;
  struct Data { Cluster_index parent; int rank; Point_index max; }; // information on a cluster
  std::vector<Data> ds_base;
  // boost::vector_property_map does resize(size+1) for every new element, don't use it
  auto ds_data=boost::make_function_property_map<std::size_t>([&ds_base](std::size_t n)->Data&{return ds_base[n];});
  auto ds_parent=boost::make_transform_value_property_map([](auto&p)->Cluster_index&{return p.parent;},ds_data);
  auto ds_rank=boost::make_transform_value_property_map([](auto&p)->int&{return p.rank;},ds_data);
  boost::disjoint_sets<decltype(ds_rank),decltype(ds_parent)> ds(ds_rank,ds_parent); // on the clusters, not directly the points
  std::vector<std::array<double,2>> persistence; // diagram (finite points)
  boost::container::flat_map<Cluster_index,Cluster_index> adj_clusters; // first: the merged cluster, second: the raw cluster
  // we only care about the raw cluster, we could use a vector to store the second, store first into a set, and only insert in the vector if merged is absent from the set

  for(Point_index i = 0; i < num_points; ++i){
    // auto&& ngb = neighbors[order[i]];
    // TODO: have a specialization where we don't need python types and py::cast
    // TODO: move py::cast and getint into Neighbors
    py::object ngb = neighbors[py::cast(order[i])]; // auto&& also seems to work
    adj_clusters.clear();
    Point_index j = i; // highest neighbor
    //for(Point_index k : ngb)
    for(auto k_py : ngb)
    {
      Point_index k = rorder[getint(k_py)];
      if(k >= i || k < 0) // ???
	continue;
      if(k < j)
	j = k;
      Cluster_index rk=raw_cluster[k];
      adj_clusters.emplace(ds.find_set(rk),rk);
      // does not insert if ck=ds.find_set(rk) already seen
      // which rk we keep from those with the same ck is arbitrary
    }
    assert((Point_index)raw_cluster.size()==i);
    if(i==j){ // local maximum -> new cluster
      Cluster_index c = n_raw_clusters++;
      ds_base.emplace_back(); // could be integrated in ds_data, but then we would check the size for every access
      ds.make_set(c);
      ds_base[c].max=i; // max
      raw_cluster.push_back(c);
      continue;
    }else{ // add i to the cluster of j
      assert(j<i);
      raw_cluster.push_back(raw_cluster[j]);
      // FIXME: we are adding point i to the raw cluster of j, but that one might not be in adj_clusters, so we may merge clusters A and B through a point of C. It is strange, but I don't know if it can really cause problems. We could just not set j at all and use arbitrarily the first element of adj_clusters.
    }
    // possibly merge clusters
    // we could sort, in case there are several merges, but that doesn't seem so useful
    Cluster_index rj = raw_cluster[j];
    Cluster_index cj = ds.find_set(rj);
    Cluster_index orig_cj = cj;
    for(auto ckk : adj_clusters){
      Cluster_index rk = ckk.second;
      Cluster_index ck = ckk.first;
      if(ck==orig_cj)
	continue;
      assert(ck==ds.find_set(rk));
      Point_index j = ds_base[cj].max;
      Point_index k = ds_base[ck].max;
      Point_index young = std::max(j,k);
      Point_index old = std::min(j,k);
      auto d_young = density[order[young]];
      auto d_i = density[order[i]];
      assert(d_young >= d_i);
      // Always merge (the non-hierarchical algorithm would only conditionally merge here
      persistence.push_back({d_young, d_i});
      assert(ds.find_set(rj) != ds.find_set(rk));
      ds.link(cj,ck);
      cj = ds.find_set(cj);
      ds_base[cj].max=old; // just one parent, no need for find_set
      // record the raw clusters, we don't know what will have already been merged.
      merges.push_back({rj, rk, d_young-d_i});
    }
  }
  {
    boost::counting_iterator<int> b(0), e(ds_base.size());
    ds.compress_sets(b, e);
    // Now we stop using find_sets and look at the parent directly
    // rank is reused to rename clusters contiguously 0, 1, etc
  }
  // Maximum for each connected component
  std::vector<double> max_cc;
  //for(Point_index i = 0; i < num_points; ++i){
  //  if(ds_base[ds_base[raw_cluster[i]].parent].max == i)
  //    max_cc.push_back(density[order[i]]);
  //}
  for(Cluster_index i = 0; i < n_raw_clusters; ++i){
    if(ds_base[i].parent == i)
      max_cc.push_back(density[order[ds_base[i].max]]);
  }
  assert((Cluster_index)(merges.size()+max_cc.size())==n_raw_clusters);

  // TODO: create a "noise" cluster, merging all those not prominent enough?
  //int nb_clusters=0;
  //for(int i=0;i<(int)ds_base.size();++i){
  //  if(ds_parent[i]!=i) continue;
  //  ds_data[i].rank=nb_clusters++;
  //}

  // Replay the merges, in increasing order of prominence, to build the hierarchy
  std::sort(merges.begin(), merges.end(), [](Merge const&a, Merge const&b){return a.persist < b.persist;});
  std::vector<std::array<Cluster_index,2>> children; children.reserve(merges.size());
  {
    struct Dat { Cluster_index parent; int rank; Cluster_index name; };
    std::vector<Dat> ds_bas(2*n_raw_clusters-1);
    Cluster_index i;
    auto ds_dat=boost::make_function_property_map<std::size_t>([&ds_bas](std::size_t n)->Dat&{return ds_bas[n];});
    auto ds_par=boost::make_transform_value_property_map([](auto&p)->Cluster_index&{return p.parent;},ds_dat);
    auto ds_ran=boost::make_transform_value_property_map([](auto&p)->int&{return p.rank;},ds_dat);
    boost::disjoint_sets<decltype(ds_ran),decltype(ds_par)> ds(ds_ran,ds_par);
    for(i=0;i<n_raw_clusters;++i) { ds.make_set(i); ds_bas[i].name=i; }
    for(Merge const& m : merges){
      Cluster_index j = ds.find_set(m.first);
      Cluster_index k = ds.find_set(m.second);
      assert(j!=k);
      children.push_back({ds_bas[j].name,ds_bas[k].name});
      ds.make_set(i);
      ds.link(i,j);
      ds.link(i,k);
      ds_bas[ds.find_set(i)].name=i;
      ++i;
    }
  }

  std::vector<Cluster_index> raw_cluster_ordered(num_points);
  for(int i=0; i<num_points; ++i) raw_cluster_ordered[i]=raw_cluster[rorder[i]];
  // return raw_cluster, children, persistence
  // TODO avoid copies: https://github.com/pybind/pybind11/issues/1042
  return py::make_tuple(
      py::array(raw_cluster_ordered.size(), raw_cluster_ordered.data()),
      py::array(children.size(), children.data()),
      py::array(persistence.size(), persistence.data()),
      py::array(max_cc.size(), max_cc.data()));
}

auto merge(py::array_t<Cluster_index, py::array::c_style> children, Cluster_index n_leaves, Cluster_index n_final){
  // Should this really be an error?
  if(n_final > n_leaves)
    throw std::runtime_error("The number of clusters required is larger than the number of mini-clusters");
  py::buffer_info cbuf = children.request();
  if((cbuf.ndim!=2 || cbuf.shape[1]!=2) && (cbuf.ndim!=1 || cbuf.shape[0]!=0))
    throw std::runtime_error("internal error: children have to be (n,2) or empty");
  const int n_merges=cbuf.shape[0];
  Cluster_index*d=(Cluster_index*)cbuf.ptr;
  // Should this really be an error?
  //std::cerr << "n_merges: " << n_merges << ", n_final: " << n_final << ", n_leaves: " << n_leaves << '\n';
  if(n_merges + n_final < n_leaves)
    throw std::runtime_error(std::string("The number of clusters required ")+std::to_string(n_final)+" is smaller than the number of connected components "+std::to_string(n_leaves-n_merges));

  struct Dat { Cluster_index parent; int rank; int name; };
  std::vector<Dat> ds_bas(2*n_leaves-1);
  auto ds_dat=boost::make_function_property_map<std::size_t>([&ds_bas](std::size_t n)->Dat&{return ds_bas[n];});
  auto ds_par=boost::make_transform_value_property_map([](auto&p)->Cluster_index&{return p.parent;},ds_dat);
  auto ds_ran=boost::make_transform_value_property_map([](auto&p)->int&{return p.rank;},ds_dat);
  boost::disjoint_sets<decltype(ds_ran),decltype(ds_par)> ds(ds_ran,ds_par);
  Cluster_index i;
  for(i=0;i<n_leaves;++i) { ds.make_set(i); ds_bas[i].name=-1; }
  for(Cluster_index m=0; m < n_leaves-n_final; ++m){
    Cluster_index j = ds.find_set(d[2*m]);
    Cluster_index k = ds.find_set(d[2*m+1]);
    assert(j!=k);
    ds.make_set(i);
    ds.link(i,j);
    ds.link(i,k);
    ++i;
  }
  Cluster_index next_cluster_name = 0;
  std::vector<Cluster_index> ret; ret.reserve(n_leaves);
  for(Cluster_index j=0;j<n_leaves;++j){
    Cluster_index k = ds.find_set(j);
    if(ds_bas[k].name == -1) ds_bas[k].name = next_cluster_name++;
    ret.push_back(ds_bas[k].name);
  }
  return py::array(ret.size(), ret.data());
}

// Do a special version when ngb is a numpy array, where we can cast to int[k][n] ?
auto plouf(py::handle ngb, py::array_t<double, py::array::c_style | py::array::forcecast> density){
  // used to be py::iterable ngb, but that's inconvenient if it doesn't come pre-sorted
  // use py::handle and check if [] (aka __getitem__) works? But then we need to build an object to pass it to []
  // (I _think_ handle is ok and we don't need object here)
  py::buffer_info wbuf = density.request();
  if(wbuf.ndim!=1)
    throw std::runtime_error("density must be 1D");
  const int n=wbuf.shape[0];
  double*d=(double*)wbuf.ptr;
  // Vector { 0, 1, ..., n-1 }
  std::vector<Point_index> order(boost::counting_iterator<Point_index>(0), boost::counting_iterator<Point_index>(n));
  // Permutation of the indices to get points in decreasing order of density
  std::sort(std::begin(order), std::end(order), [=](Point_index i, Point_index j){ return d[i] > d[j]; });
  // Inverse permutation
  std::vector<Point_index> rorder(n);
  for(Point_index i : boost::irange(0, n)) rorder[order[i]] = i;
  // Used as:
  // order[i] is the index of the point with i-th largest density
  // rorder[i] is the rank of the i-th point in order of decreasing density
  // TODO: put a wrapper on ngb and d so we don't need to pass (r)order (there is still the issue of reordering the output)
  return tomato(n, ngb, d, order, rorder);
}
#if 0
struct Vec {
  int const* base;
  int k;
  int const* begin()const{return base;}
  int const* end()const{return base+k;}
  Vec operator*()const{return *this;}
};
auto plaf(py::array_t<int, py::array::c_style | py::array::forcecast> ngb, py::array_t<double, py::array::c_style | py::array::forcecast> density, double threshold){
  py::buffer_info nbuf = ngb.request();
  if(nbuf.ndim!=2)
    throw std::runtime_error("neighbors must be 2D");
  const int n=nbuf.shape[0];
  const int k=nbuf.shape[1];
  int*nptr=(int*)nbuf.ptr;
  auto neighbors=boost::adaptors::transform(boost::irange(0,n),[=](int i){return Vec{nptr+i*k,k};});
  py::buffer_info wbuf = density.request();
  if(wbuf.ndim!=1)
    throw std::runtime_error("density must be 1D");
  if(n!=wbuf.shape[0])
    throw std::runtime_error("incompatible sizes");
  double*d=(double*)wbuf.ptr;
  return tomato(n,neighbors,d);
}
#endif
PYBIND11_MODULE(_tomato, m) {
  m.doc() = "Internals of tomato clustering";
  m.def("doit", &plouf, "does the clustering");
  //m.def("doit2", &plaf, "does the clustering faster");
  m.def("merge", &merge, "merge clusters");
}

// https://github.com/pybind/pybind11/issues/1042 pour convertir vector en numpy array
//
// py::isinstance<py::array_t<std::int32_t>> (ou py::isinstance<py::array> et tester dtype) et flags&c_style
// ou overload (en virant forcecast?)
// aussi le faire au cas où on n'aurait pas un tableau, mais où chaque liste de voisins serait un tableau ?
