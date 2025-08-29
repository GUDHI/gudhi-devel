/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - 2025/03 Vincent Rouvreau & Hannah Schreiber: Use nanobind instead of PyBind11 for python bindings.
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>
#include <iostream>  //std::cerr
#include <array>

#include <boost/container/flat_map.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/transform_value_property_map.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <python_interfaces/numpy_utils.h>

namespace nb = nanobind;

template <class T, class = std::enable_if_t<std::is_integral<T>::value>>
int getint(int i)
{
  return i;
}

// Gcc-8 has a bug that breaks this version, fixed in gcc-9
// template<class T, class=decltype(std::declval<T>().template cast<int>())>
// int getint(T i){return i.template cast<int>();}
template <class T>
auto getint(T i) -> decltype(nb::cast<int>(i))
{
  return nb::cast<int>(i);
}

// Raw clusters are clusters obtained through single-linkage, no merging.

typedef int Point_index;
typedef int Cluster_index;

struct Merge {
  Cluster_index first, second;
  double persist;
};

template <class Density_view>
auto tomato(Point_index num_points,
            nb::object const& neighbors,
            const Density_view& density,
            std::vector<Point_index> const& order,
            std::vector<Point_index> const& rorder)
{
  // point index --> index of raw cluster it belongs to
  std::vector<Cluster_index> raw_cluster;
  raw_cluster.reserve(num_points);
  // cluster index --> index of top point in the cluster
  Cluster_index n_raw_clusters = 0;  // current number of raw clusters seen
  //
  std::vector<Merge> merges;

  struct Data {
    Cluster_index parent;
    int rank;
    Point_index max;
  };  // information on a cluster

  std::vector<Data> ds_base;
  // boost::vector_property_map does resize(size+1) for every new element, don't use it
  auto ds_data =
      boost::make_function_property_map<std::size_t>([&ds_base](std::size_t n) -> Data& { return ds_base[n]; });
  auto ds_parent =
      boost::make_transform_value_property_map([](auto& p) -> Cluster_index& { return p.parent; }, ds_data);
  auto ds_rank = boost::make_transform_value_property_map([](auto& p) -> int& { return p.rank; }, ds_data);
  // on the clusters, not directly the points
  boost::disjoint_sets<decltype(ds_rank), decltype(ds_parent)> ds(ds_rank, ds_parent);
  // diagram (finite points)
  std::vector<std::array<double, 2> > persistence;
  // first: the merged cluster, second: the raw cluster
  // we only care about the raw cluster, we could use a vector to store the second, store first into a set, and only
  // insert in the vector if merged is absent from the set
  boost::container::flat_map<Cluster_index, Cluster_index> adj_clusters;

  for (Point_index i = 0; i < num_points; ++i) {
    // auto&& ngb = neighbors[order[i]];
    // TODO: have a specialization where we don't need python types and py::cast
    // TODO: move py::cast and getint into Neighbors
    nb::object ngb = neighbors[order[i]];  // auto&& also seems to work
    adj_clusters.clear();
    Point_index j = i;  // highest neighbor
    // for(Point_index k : ngb)
    for (auto k_py : ngb) {
      Point_index k = rorder[getint(k_py)];
      if (k >= i || k < 0)  // ???
        continue;
      if (k < j) j = k;
      Cluster_index rk = raw_cluster[k];
      adj_clusters.emplace(ds.find_set(rk), rk);
      // does not insert if ck=ds.find_set(rk) already seen
      // which rk we keep from those with the same ck is arbitrary
    }
    assert((Point_index)raw_cluster.size() == i);
    if (i == j) {  // local maximum -> new cluster
      Cluster_index c = n_raw_clusters++;
      ds_base.emplace_back();  // could be integrated in ds_data, but then we would check the size for every access
      ds.make_set(c);
      ds_base[c].max = i;  // max
      raw_cluster.push_back(c);
      continue;
    } else {  // add i to the cluster of j
      assert(j < i);
      raw_cluster.push_back(raw_cluster[j]);
      // FIXME: we are adding point i to the raw cluster of j, but that one might not be in adj_clusters, so we may
      // merge clusters A and B through a point of C. It is strange, but I don't know if it can really cause problems.
      // We could just not set j at all and use arbitrarily the first element of adj_clusters.
    }
    // possibly merge clusters
    // we could sort, in case there are several merges, but that doesn't seem so useful
    Cluster_index rj = raw_cluster[j];
    Cluster_index cj = ds.find_set(rj);
    Cluster_index orig_cj = cj;
    for (auto ckk : adj_clusters) {
      Cluster_index rk = ckk.second;
      Cluster_index ck = ckk.first;
      if (ck == orig_cj) continue;
      assert(ck == ds.find_set(rk));
      Point_index j = ds_base[cj].max;
      Point_index k = ds_base[ck].max;
      Point_index young = std::max(j, k);
      Point_index old = std::min(j, k);
      auto d_young = density(order[young]);
      auto d_i = density(order[i]);
      assert(d_young >= d_i);
      // Always merge (the non-hierarchical algorithm would only conditionally merge here
      persistence.emplace_back(std::array<double, 2>{d_young, d_i});
      assert(ds.find_set(rj) != ds.find_set(rk));
      ds.link(cj, ck);
      cj = ds.find_set(cj);
      ds_base[cj].max = old;  // just one parent, no need for find_set
      // record the raw clusters, we don't know what will have already been merged.
      merges.push_back({rj, rk, d_young - d_i});
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
  for (Cluster_index i = 0; i < n_raw_clusters; ++i) {
    if (ds_base[i].parent == i) max_cc.push_back(density(order[ds_base[i].max]));
  }
  assert((Cluster_index)(merges.size() + max_cc.size()) == n_raw_clusters);

  // TODO: create a "noise" cluster, merging all those not prominent enough?

  // Replay the merges, in increasing order of prominence, to build the hierarchy
  std::sort(merges.begin(), merges.end(), [](Merge const& a, Merge const& b) { return a.persist < b.persist; });
  std::vector<std::array<Cluster_index, 2> > children;
  children.reserve(merges.size() * 2);
  {
    struct Dat {
      Cluster_index parent;
      int rank;
      Cluster_index name;
    };

    std::vector<Dat> ds_bas(2 * n_raw_clusters - 1);
    Cluster_index i;
    auto ds_dat =
        boost::make_function_property_map<std::size_t>([&ds_bas](std::size_t n) -> Dat& { return ds_bas[n]; });
    auto ds_par = boost::make_transform_value_property_map([](auto& p) -> Cluster_index& { return p.parent; }, ds_dat);
    auto ds_ran = boost::make_transform_value_property_map([](auto& p) -> int& { return p.rank; }, ds_dat);
    boost::disjoint_sets<decltype(ds_ran), decltype(ds_par)> ds(ds_ran, ds_par);
    for (i = 0; i < n_raw_clusters; ++i) {
      ds.make_set(i);
      ds_bas[i].name = i;
    }
    for (Merge const& m : merges) {
      Cluster_index j = ds.find_set(m.first);
      Cluster_index k = ds.find_set(m.second);
      assert(j != k);
      children.emplace_back(std::array<Cluster_index, 2>{ds_bas[j].name, ds_bas[k].name});
      ds.make_set(i);
      ds.link(i, j);
      ds.link(ds.find_set(i), k);
      ds_bas[ds.find_set(i)].name = i;
      ++i;
    }
  }

  std::vector<Cluster_index> raw_cluster_ordered(num_points);
  for (int i = 0; i < num_points; ++i) raw_cluster_ordered[i] = raw_cluster[rorder[i]];

  return nb::make_tuple(_wrap_as_numpy_array(std::move(raw_cluster_ordered), raw_cluster_ordered.size()),
                        _wrap_as_numpy_array(std::move(children)),
                        _wrap_as_numpy_array(std::move(persistence)),
                        _wrap_as_numpy_array(std::move(max_cc), max_cc.size()));
}

auto merge(nb::ndarray<const Cluster_index, nb::ndim<2> > children,
           Cluster_index n_leaves,
           Cluster_index n_final)
{
  if (n_final > n_leaves) {
    std::cerr << "The number of clusters required " << n_final << " is larger than the number of mini-clusters "
              << n_leaves << '\n';
    n_final = n_leaves;  // or return something special and let Tomato use leaf_labels_?
  }

  auto children_view = children.view();

  if (children_view.shape(1) != 2) throw std::runtime_error("internal error: `children` has to be (n,2) shaped.");
  const int n_merges = children_view.shape(0);
  if (n_merges + n_final < n_leaves) {
    std::cerr << "The number of clusters required " << n_final << " is smaller than the number of connected components "
              << n_leaves - n_merges << '\n';
    n_final = n_leaves - n_merges;
  }

  struct Dat {
    Cluster_index parent;
    int rank;
    int name;
  };

  std::vector<Dat> ds_bas(2 * n_leaves - 1);
  auto ds_dat = boost::make_function_property_map<std::size_t>([&ds_bas](std::size_t n) -> Dat& { return ds_bas[n]; });
  auto ds_par = boost::make_transform_value_property_map([](auto& p) -> Cluster_index& { return p.parent; }, ds_dat);
  auto ds_ran = boost::make_transform_value_property_map([](auto& p) -> int& { return p.rank; }, ds_dat);
  boost::disjoint_sets<decltype(ds_ran), decltype(ds_par)> ds(ds_ran, ds_par);
  Cluster_index i;
  for (i = 0; i < n_leaves; ++i) {
    ds.make_set(i);
    ds_bas[i].name = -1;
  }
  for (Cluster_index m = 0; m < n_leaves - n_final; ++m) {
    Cluster_index j = ds.find_set(children_view(m, 0));
    Cluster_index k = ds.find_set(children_view(m, 1));
    assert(j != k);
    ds.make_set(i);
    ds.link(i, j);
    ds.link(ds.find_set(i), k);
    ++i;
  }
  Cluster_index next_cluster_name = 0;
  std::vector<Cluster_index> ret;
  ret.reserve(n_leaves);
  for (Cluster_index j = 0; j < n_leaves; ++j) {
    Cluster_index k = ds.find_set(j);
    if (ds_bas[k].name == -1) ds_bas[k].name = next_cluster_name++;
    ret.push_back(ds_bas[k].name);
  }
  return _wrap_as_numpy_array(std::move(ret), ret.size());
}

// TODO: Do a special version when ngb is a numpy array, where we can cast to int[k][n] ?
// also do this in the case where we don't have an array, but each list of neighbours is an array ?
auto hierarchy(nb::object ngb, nb::ndarray<const double, nb::ndim<1> > density)
{
  // used to be py::iterable ngb, but that's inconvenient if it doesn't come pre-sorted
  // use py::handle and check if [] (aka __getitem__) works? But then we need to build an object to pass it to []
  // (I _think_ handle is ok and we don't need object here)
  const int n = density.shape(0);
  // Vector { 0, 1, ..., n-1 }
  std::vector<Point_index> order(boost::counting_iterator<Point_index>(0), boost::counting_iterator<Point_index>(n));
  auto density_view = density.view();
  // Permutation of the indices to get points in decreasing order of density
  std::sort(std::begin(order), std::end(order), [&density_view](Point_index i, Point_index j) {
    return density_view(i) > density_view(j);
  });
  // Inverse permutation
  std::vector<Point_index> rorder(n);
  for (Point_index i : boost::irange(0, n)) rorder[order[i]] = i;
  // Used as:
  // order[i] is the index of the point with i-th largest density
  // rorder[i] is the rank of the i-th point in order of decreasing density
  // TODO: put a wrapper on ngb and density so we don't need to pass (r)order (there is still the issue of reordering
  // the output)
  return tomato(n, ngb, density_view, order, rorder);
}

NB_MODULE(_tomato_ext, m)
{
  m.attr("__license__") = "MIT";
  m.doc() = "Internals of tomato clustering";
  m.def("_hierarchy", &hierarchy);  // does the clustering
  m.def("_merge", &merge);          // merge clusters
}
