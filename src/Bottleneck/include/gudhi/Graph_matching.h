/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Francois Godi
 *
 *    Copyright (C) 2015  INRIA Sophia-Antipolis (France)
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

#ifndef SRC_BOTTLENECK_INCLUDE_GUDHI_GRAPH_MATCHING_H_
#define SRC_BOTTLENECK_INCLUDE_GUDHI_GRAPH_MATCHING_H_
//#define DEBUG
#include <deque>
#include <list>
#include <vector>

#include "Layered_neighbors_finder.h"

namespace Gudhi {

namespace bottleneck {

template<typename Persistence_diagram1, typename Persistence_diagram2>
double bottleneck_distance(const Persistence_diagram1& diag1, const Persistence_diagram2& diag2, double e = 0.);

class Graph_matching {
public:
    explicit Graph_matching(const Persistence_diagrams_graph& g);
    Graph_matching& operator=(const Graph_matching& m);
    bool perfect() const;
    bool multi_augment();
    void set_r(double r);

private:
    const Persistence_diagrams_graph& g;
    double r;
    std::vector<int> v_to_u;
    std::list<int> unmatched_in_u;

    std::unique_ptr<Layered_neighbors_finder> layering() const;
    bool augment(Layered_neighbors_finder & layered_nf, int u_start_index, int max_depth);
    void update(std::deque<int> & path);
};

Graph_matching::Graph_matching(const Persistence_diagrams_graph& g)
    : g(g), r(0), v_to_u(g.size(), null_point_index()), unmatched_in_u() {
    for (int u_point_index = 0; u_point_index < g.size(); ++u_point_index)
        unmatched_in_u.emplace_back(u_point_index);
}

Graph_matching& Graph_matching::operator=(const Graph_matching& m) {
    r = m.r;
    v_to_u = m.v_to_u;
    unmatched_in_u = m.unmatched_in_u;
    return *this;
}

/* inline */ bool Graph_matching::perfect() const {
#ifdef DEBUG
    std::cout << "  perfect? unmatched_in_u.size = " << unmatched_in_u.size() << std::endl << std::flush;
#endif
    return unmatched_in_u.empty();
}

/* inline */ bool Graph_matching::multi_augment() {
    if (perfect())
        return false;
#ifdef DEBUG
    std::cout << "  multi augment" << std::endl << std::flush;
#endif
    Layered_neighbors_finder layered_nf = *layering();
    double rn = sqrt(g.size());
    int max_depth = layered_nf.vlayers_number()*2 - 1;
#ifdef DEBUG
    std::cout<< "    nb_max_layer = " << max_depth << std::endl << std::flush;
#endif
    // verification of a necessary criterion
    if (max_depth <0 || (unmatched_in_u.size() > rn && max_depth >= rn))
        return false;
    bool successful = false;
    std::list<int> tries(unmatched_in_u);
    for (auto it = tries.cbegin(); it != tries.cend(); it++){
        const bool tmp = augment(layered_nf, *it, max_depth); //force augment evaluation even if successful is already true
        successful = successful || tmp;
    }
    return successful;
}

/* inline */ void Graph_matching::set_r(double r) {
    this->r = r;
}

bool Graph_matching::augment(Layered_neighbors_finder & layered_nf, int u_start_index, int max_depth) {
#ifdef DEBUG
    std::cout << "    augment" << std::endl << std::flush;
#endif
#ifdef DEBUG
    std::cout << "      u_start_index = " << u_start_index << std::endl << std::flush;
#endif
    std::deque<int> path;
    path.emplace_back(u_start_index);
    // u_start is a point from U
    do {
#ifdef DEBUG
        std::cout << "      do" << std::endl << std::flush;
        std::cout << "        path.size = " << static_cast<int>(path.size()) << std::endl << std::flush;
#endif
        if (static_cast<int>(path.size()) > max_depth) {
            path.pop_back();
            path.pop_back();
        }
        if (path.empty())
            return false;
        path.emplace_back(layered_nf.pull_near(path.back(), static_cast<int>(path.size())/2));
        while (path.back() == null_point_index()) {
            path.pop_back();
            path.pop_back();
            if (path.empty())
                return false;
            path.pop_back();
            path.emplace_back(layered_nf.pull_near(path.back(), path.size() / 2));
        }
        path.emplace_back(v_to_u.at(path.back()));
#ifdef DEBUG
        std::cout << "        v_to_u = " << path.back() << std::endl << std::flush;
#endif
    } while (path.back() != null_point_index());
    path.pop_back();
    update(path);
    return true;
}

std::unique_ptr<Layered_neighbors_finder> Graph_matching::layering() const {
#ifdef DEBUG
    std::cout << "    layering" << std::endl << std::flush;
#endif
    bool end = false;
    int layer = 0;
    std::list<int> u_vertices(unmatched_in_u);
    std::list<int> v_vertices;
    Neighbors_finder nf(g, r);
    for (int v_point_index = 0; v_point_index < g.size(); ++v_point_index)
        nf.add(v_point_index);
    std::unique_ptr<Layered_neighbors_finder> layered_nf(new Layered_neighbors_finder(g, r));
    while (!u_vertices.empty()) {
        for (auto it = u_vertices.cbegin(); it != u_vertices.cend(); ++it) {
            std::unique_ptr< std::list<int> > u_succ = std::move(nf.pull_all_near(*it));
            for (auto it = u_succ->cbegin(); it != u_succ->cend(); ++it) {
                layered_nf->add(*it, layer);
                v_vertices.emplace_back(*it);
            }
        }
        u_vertices.clear();
        for (auto it = v_vertices.cbegin(); it != v_vertices.cend(); it++) {
            if (v_to_u.at(*it) == null_point_index())
                end = true;
            else
                u_vertices.emplace_back(v_to_u.at(*it));
        }
        if (end)
            return layered_nf;
        v_vertices.clear();
        layer++;
    }
    return layered_nf;
}

void Graph_matching::update(std::deque<int>& path) {
#ifdef DEBUG
    std::cout << "      update" << std::endl << std::flush;
#endif
    unmatched_in_u.remove(path.front());
    for (auto it = path.cbegin(); it != path.cend(); ++it) {
        int tmp = *it;
        ++it;
        v_to_u[*it] = tmp;
    }
}

template<typename Persistence_diagram1, typename Persistence_diagram2>
double bottleneck_distance(const Persistence_diagram1 &diag1, const Persistence_diagram2 &diag2, double e) {
    Persistence_diagrams_graph g(diag1, diag2, e);
    std::unique_ptr< std::vector<double> > sd = std::move(g.sorted_distances());
    int idmin = 0;
    int idmax = sd->size() - 1;
    double alpha = pow(sd->size(), 0.25);
    Graph_matching m(g);
    Graph_matching biggest_unperfect(g);
    while (idmin != idmax) {
        int pas = static_cast<int>((idmax - idmin) / alpha);
        m.set_r(sd->at(idmin + pas));
        while (m.multi_augment());
        if (m.perfect()) {
            idmax = idmin + pas;
            m = biggest_unperfect;
        } else {
            biggest_unperfect = m;
            idmin = idmin + pas + 1;
        }
    }
    double b = sd->at(idmin);
    return b;
}

}  // namespace bottleneck

}  // namespace Gudhi

#endif  // SRC_BOTTLENECK_INCLUDE_GUDHI_GRAPH_MATCHING_H_
