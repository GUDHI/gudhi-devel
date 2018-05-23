#ifndef LAZY_TOPLEX_MAP_H
#define LAZY_TOPLEX_MAP_H

#include <gudhi/Toplex_map.h>
#include <boost/heap/fibonacci_heap.hpp>
#include <cmath>

namespace Gudhi {

class Lazy_Toplex_map {

public:

    /** Vertex is the type of vertices.
     * \ingroup toplex_map   */
    typedef Toplex_map::Vertex Vertex;

    /** Simplex is the type of simplices.
     * \ingroup toplex_map   */
    typedef Toplex_map::Simplex Simplex;

    /** The type of the pointers to maximal simplices.
     * \ingroup toplex_map   */
    typedef Toplex_map::Simplex_ptr Simplex_ptr;

    /** The type of the sets of Simplex_ptr.
     * \ingroup toplex_map   */
    typedef Toplex_map::Simplex_ptr_set Simplex_ptr_set;

    template <typename Input_vertex_range>
    void insert_max_simplex(const Input_vertex_range &vertex_range);
    template <typename Input_vertex_range>
    bool insert_simplex(const Input_vertex_range &vertex_range);
    template <typename Input_vertex_range>
    void remove_simplex(const Input_vertex_range &vertex_range);

    template <typename Input_vertex_range>
    bool membership(const Input_vertex_range &vertex_range);
    template <typename Input_vertex_range>
    bool all_facets_inside(const Input_vertex_range &vertex_range);

    Vertex contraction(const Vertex x, const Vertex y);

    std::size_t num_simplices() const;

    std::unordered_map<Vertex, std::size_t> gamma0_lbounds;

private:
    template <typename Input_vertex_range>
    void erase_max(const Input_vertex_range &vertex_range);
    template <typename Input_vertex_range>
    Vertex best_index(const Input_vertex_range &vertex_range);
    void clean(const Vertex v);

    std::unordered_map<Vertex, Simplex_ptr_set> t0;
    bool empty_toplex; // Is the empty simplex a toplex ?

    typedef boost::heap::fibonacci_heap<std::pair<std::size_t,Vertex>> PriorityQueue;
    PriorityQueue cleaning_priority;
    std::unordered_map<Vertex, PriorityQueue::handle_type> cp_handles;

    std::size_t get_gamma0_lbound(const Vertex v) const;

    std::size_t size_lbound = 0;
    std::size_t size = 0;

    const double alpha = 4; //time
    const double betta = 8; //memory
};

template <typename Input_vertex_range>
void Lazy_Toplex_map::insert_max_simplex(const Input_vertex_range &vertex_range){
    for(const Vertex& v : vertex_range)
        if(!gamma0_lbounds.count(v)) gamma0_lbounds.emplace(v,1);
        else gamma0_lbounds[v]++;
    size_lbound++;
    insert_simplex(vertex_range);
}

template <typename Input_vertex_range>
bool Lazy_Toplex_map::insert_simplex(const Input_vertex_range &vertex_range){
    Simplex sigma(vertex_range.begin(),vertex_range.end());
    empty_toplex = (sigma.size()==0); //vérifier la gestion de empty face
    Simplex_ptr sptr = std::make_shared<Simplex>(sigma);
    bool inserted = false;
    for(const Vertex& v : sigma){
        if(!t0.count(v)){
            t0.emplace(v, Simplex_ptr_set());
            auto v_handle = cleaning_priority.push(std::make_pair(0, v));
            cp_handles.emplace(v, v_handle);
        }
        inserted = t0.at(v).emplace(sptr).second;
        cleaning_priority.update(cp_handles.at(v), std::make_pair(t0.at(v).size() - get_gamma0_lbound(v),v));
    }
    if(inserted)
        size++;
    if(size > (size_lbound+1) * betta)
        clean(cleaning_priority.top().second);
    return inserted;
}

template <typename Input_vertex_range>
void Lazy_Toplex_map::remove_simplex(const Input_vertex_range &vertex_range){
    if(vertex_range.begin()==vertex_range.end()){
        t0.clear();
        gamma0_lbounds.clear();
        cleaning_priority.clear();
        size_lbound = 0;
        size = 0;
        empty_toplex = false;
    }
    else {
        const Vertex& v = best_index(vertex_range);
        //Copy constructor needed because the set is modified
        if(t0.count(v)) for(const Simplex_ptr& sptr : Simplex_ptr_set(t0.at(v)))
            if(included(vertex_range, *sptr)){
                erase_max(*sptr);
                for(const Simplex& f : facets(vertex_range))
                    insert_max_simplex(f);
            }
    }
}

template <typename Input_vertex_range>
bool Lazy_Toplex_map::membership(const Input_vertex_range &vertex_range){
    if(t0.size()==0 && !empty_toplex) return false; //empty complex
    if(vertex_range.begin()==vertex_range.end()) return true; //empty query simplex
    Vertex v = best_index(vertex_range);
    if(!t0.count(v))  return false;
    for(const Simplex_ptr& sptr : t0.at(v))
        if(included(vertex_range, *sptr)) return true;
    return false;
}

template <typename Input_vertex_range>
bool Lazy_Toplex_map::all_facets_inside(const Input_vertex_range &vertex_range){
    Simplex sigma(vertex_range.begin(),vertex_range.end());
    Vertex v = best_index(sigma);
    if(!t0.count(v))  return false;
    Simplex f = sigma; f.erase(v);
    if(!membership(f)) return false;
    std::unordered_set<Vertex> facets_inside;
    for(const Simplex_ptr& sptr : t0.at(v))
        for(const Vertex& w : sigma){
            f = sigma; f.erase(w);
            if(included(f, *sptr)) facets_inside.insert(w);
        }
    return facets_inside.size() == sigma.size() - 1;
}

/* Returns the remaining vertex */
Toplex_map::Vertex Lazy_Toplex_map::contraction(const Vertex x, const Vertex y){
    if(!t0.count(x)) return y;
    if(!t0.count(y)) return x;
    Vertex k, d;
    if(t0.at(x).size() > t0.at(y).size())
        k=x, d=y;
    else
        k=y, d=x;
    //Copy constructor needed because the set is modified
    for(const Simplex_ptr& sptr : Simplex_ptr_set(t0.at(d))){
        Simplex sigma(*sptr);
        erase_max(sigma);
        sigma.erase(d);
        sigma.insert(k);
        insert_simplex(sigma);
    }
    t0.erase(d);
    return k;
}

/* No facets insert_simplexed */
template <typename Input_vertex_range>
inline void Lazy_Toplex_map::erase_max(const Input_vertex_range &vertex_range){
    Simplex sigma(vertex_range.begin(),vertex_range.end());
    empty_toplex = false;
    Simplex_ptr sptr = std::make_shared<Simplex>(sigma);
    bool erased=false;
    for(const Vertex& v : sigma){
        erased = t0.at(v).erase(sptr) > 0;
        if(t0.at(v).size()==0)
            t0.erase(v);
    }
    if (erased)
        size--;
}

template <typename Input_vertex_range>
Toplex_map::Vertex Lazy_Toplex_map::best_index(const Input_vertex_range &vertex_range){
    Simplex tau(vertex_range.begin(),vertex_range.end());
    std::size_t min = std::numeric_limits<size_t>::max(); Vertex arg_min = -1;
    for(const Vertex& v : tau)
        if(!t0.count(v)) return v;
        else if(t0.at(v).size() < min)
            min = t0.at(v).size(), arg_min = v;
    if(min > alpha * get_gamma0_lbound(arg_min))
        clean(arg_min);
    return arg_min;
}

std::size_t Lazy_Toplex_map::get_gamma0_lbound(const Vertex v) const{
    return gamma0_lbounds.count(v) ? gamma0_lbounds.at(v) : 0;
}


void Lazy_Toplex_map::clean(const Vertex v){
    Toplex_map toplices;
    std::unordered_map<int, std::vector<Simplex>> dsorted_simplices;
    std::size_t max_dim = 0;
    for(const Simplex_ptr& sptr : Simplex_ptr_set(t0.at(v))){
        if(sptr->size() > max_dim){
            for(std::size_t d = max_dim+1; d<=sptr->size(); d++)
                dsorted_simplices.emplace(d, std::vector<Simplex>());
            max_dim = sptr->size();
        }
        dsorted_simplices[sptr->size()].emplace_back(*sptr);
        erase_max(*sptr);
    }
    for(std::size_t d = max_dim; d>=1; d--)
        for(const Simplex &s : dsorted_simplices.at(d))
            if(!toplices.membership(s))
                toplices.insert_independent_simplex(s);
    Simplex sv; sv.insert(v);
    auto clean_cofaces = toplices.maximal_cofaces(sv);
    size_lbound = size_lbound - get_gamma0_lbound(v) + clean_cofaces.size();
    gamma0_lbounds[v] = clean_cofaces.size();
    for(const Simplex_ptr& sptr : clean_cofaces)
        insert_simplex(*sptr);
}

std::size_t Lazy_Toplex_map::num_simplices() const{
    return size;
}

} //namespace Gudhi

#endif /* LAZY_TOPLEX_MAP_H */
