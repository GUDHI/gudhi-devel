#ifndef TOPLEX_MAP_H
#define TOPLEX_MAP_H

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>

#define vertex_upper_bound std::numeric_limits<Vertex>::max()

namespace Gudhi {

/** Vertex is the type of vertices.
 * \ingroup toplex_map   */
typedef std::size_t Vertex;

/** Simplex is the type of simplices.
 * \ingroup toplex_map   */
typedef std::unordered_set<Vertex> Simplex;

/** A Toplex_map represents the simplicial complex.
 * A "toplex" is a maximal simplex.
 * \ingroup toplex_map   */
class Toplex_map {

public:
    /** The type of the pointers to maximal simplices.
     * \ingroup toplex_map   */
    typedef std::shared_ptr<Simplex> Simplex_ptr;

    struct Sptr_hash{ std::size_t operator()(const Simplex_ptr& s) const; };
    struct Sptr_equal{ std::size_t operator()(const Simplex_ptr& a, const Simplex_ptr& b) const; };
    /** The type of the sets of Simplex_ptr.
     * \ingroup toplex_map   */
    typedef std::unordered_set<Simplex_ptr, Sptr_hash, Sptr_equal> Simplex_ptr_set;

    /** \brief Adds the given simplex to the complex.
     * Nothing happens if the simplex has a coface in the complex.
     * \ingroup toplex_map   */
    template <typename Input_vertex_range>
    void insert_simplex(const Input_vertex_range &vertex_range);

    /** \brief Removes the given simplex and its cofaces from the complex.
     * Its faces are kept inside.
     * \ingroup toplex_map   */
    template <typename Input_vertex_range>
    void remove_simplex(const Input_vertex_range &vertex_range);

    /** Does a simplex belong to the complex ?
     * \ingroup toplex_map   */
    template <typename Input_vertex_range>
    bool membership(const Input_vertex_range &vertex_range) const;

    /** Does a simplex is a toplex ?
     * \ingroup toplex_map   */
    template <typename Input_vertex_range>
    bool maximality(const Input_vertex_range &vertex_range) const;

    /** Gives a set of pointers to the maximal cofaces of a simplex.
     * Gives the toplices if given the empty simplex.
     * Gives not more than max_number maximal cofaces if max_number is strictly positive.
     * \ingroup toplex_map   */
    template <typename Input_vertex_range>
    Simplex_ptr_set maximal_cofaces(const Input_vertex_range &vertex_range, const std::size_t max_number = 0) const;

    /** Contracts one edge in the complex.
     * The edge has to verify the link condition if you want to preserve topology.
     * Returns the remaining vertex.
     * \ingroup toplex_map   */
    Vertex contraction(const Vertex x, const Vertex y);

    /** Adds the given simplex to the complex.
     * The simplex must not have neither maximal face nor coface in the complex.
     * \ingroup toplex_map   */
    template <typename Input_vertex_range>
    void insert_independent_simplex(const Input_vertex_range &vertex_range);


    /** \internal Removes a toplex without adding facets after.
     * \ingroup toplex_map   */
    void erase_maximal(const Simplex_ptr& sptr);

    /** Removes a vertex from any simplex containing it.
     * \ingroup toplex_map   */
    void remove_vertex(const Vertex x);

    /** \brief Number of maximal simplices.
         * /!\ Not efficient !
         * \ingroup toplex_map   */
    std::size_t num_simplices() const;

protected:
    /** \internal Gives an index in order to look for a simplex quickly.
     * \ingroup toplex_map   */
    template <typename Input_vertex_range>
    Vertex best_index(const Input_vertex_range &vertex_range) const;
    
    /** \internal The map from vertices to toplices
     * \ingroup toplex_map   */
    std::unordered_map<Vertex, Simplex_ptr_set> t0;

};

typedef Toplex_map::Simplex_ptr Simplex_ptr;
typedef Toplex_map::Simplex_ptr_set Simplex_ptr_set;

// Pointers are also used as key in the hash sets.
template <typename Input_vertex_range>
Simplex_ptr get_key(const Input_vertex_range &vertex_range);

// Is the first simplex a face of the second ?
template <typename Input_vertex_range1, typename Input_vertex_range2>
bool included(const Input_vertex_range1 &vertex_range1, const Input_vertex_range2 &vertex_range2);

// All the facets of the given simplex.
template <typename Input_vertex_range>
std::vector<Simplex> facets(const Input_vertex_range &vertex_range);

template <typename Input_vertex_range>
void Toplex_map::insert_simplex(const Input_vertex_range &vertex_range){
    if(membership(vertex_range)) return;
    bool replace_facets = true;
    for(const Simplex& facet : facets(vertex_range))
        if(!maximality(facet))
        {
            replace_facets=false;
            break;
        }
    if(replace_facets)
        for(const Simplex& facet : facets(vertex_range))
            erase_maximal(get_key(facet));
    else
        for(const Vertex& v : vertex_range)
            if(t0.count(v))  for(const Simplex_ptr& fptr : Simplex_ptr_set(t0.at(v)))
                //Copy constructor needed because the set is modified
                if(included(*fptr,vertex_range)) erase_maximal(fptr);
    // We erase all the maximal faces of the simplex
    insert_independent_simplex(vertex_range);
}

template <typename Input_vertex_range>
void Toplex_map::remove_simplex(const Input_vertex_range &vertex_range){
    if(vertex_range.begin()==vertex_range.end())
        t0.clear();
    // Removal of the empty simplex means cleaning everything
    else {
        const Vertex& v = best_index(vertex_range);
        if(t0.count(v)) for(const Simplex_ptr& sptr : Simplex_ptr_set(t0.at(v)))
            //Copy constructor needed because the set is modified
            if(included(vertex_range, *sptr)){
                erase_maximal(sptr);
                for(const Simplex& f : facets(vertex_range))
                    if(!membership(f)) insert_independent_simplex(f);
                // We add the facets which are new maximal simplices
            }
    }
}

template <typename Input_vertex_range>
bool Toplex_map::membership(const Input_vertex_range &vertex_range) const{
    if(t0.size()==0) return false;
    const Vertex& v = best_index(vertex_range);
    if(!t0.count(v))  return false;
    if(maximality(vertex_range)) return true;
    for(const Simplex_ptr& sptr : t0.at(v))
        if(included(vertex_range, *sptr))
            return true;
    return false;
}

template <typename Input_vertex_range>
bool Toplex_map::maximality(const Input_vertex_range &vertex_range) const{
    const Vertex& v =  best_index(vertex_range);
    if(!t0.count(v)) return false;
    return t0.at(v).count(get_key(vertex_range));
}

template <typename Input_vertex_range>
Simplex_ptr_set Toplex_map::maximal_cofaces(const Input_vertex_range &vertex_range, const std::size_t max_number) const{
    Simplex_ptr_set cofaces;
    if(maximality(vertex_range))
        cofaces.emplace(get_key(vertex_range));
    else if(vertex_range.begin()==vertex_range.end())
        for(const auto& kv : t0)
            for(const Simplex_ptr& sptr : kv.second){
                //kv.second is a Simplex_ptr_set
                cofaces.emplace(sptr);
                if(cofaces.size()==max_number)
                    return cofaces;
            }
    else {
        const Vertex& v = best_index(vertex_range);
        if(t0.count(v)) for(const Simplex_ptr& sptr : t0.at(v))
            if(included(vertex_range, *sptr)){
                cofaces.emplace(sptr);
                if(cofaces.size()==max_number)
                    return cofaces;
            }
    }
    return cofaces;
}

Vertex Toplex_map::contraction(const Vertex x, const Vertex y){
    if(!t0.count(x)) return y;
    if(!t0.count(y)) return x;
    int k, d;
    if(t0.at(x).size() > t0.at(y).size())
        k=x, d=y;
    else
        k=y, d=x;
    for(const Simplex_ptr& sptr : Simplex_ptr_set(t0.at(d))){
        //Copy constructor needed because the set is modified
        Simplex sigma(*sptr);
        erase_maximal(sptr);
        sigma.erase(d);
        sigma.insert(k);
        insert_simplex(sigma);
    }
    return k;
}

template <typename Input_vertex_range>
void Toplex_map::insert_independent_simplex(const Input_vertex_range &vertex_range){
    for(const Vertex& v : vertex_range){
        if(!t0.count(v)) t0.emplace(v, Simplex_ptr_set());
        t0.at(v).emplace(get_key(vertex_range));
    }
}

void Toplex_map::remove_vertex(const Vertex x){
    for(const Simplex_ptr& sptr : Simplex_ptr_set(t0.at(x))){
        Simplex sigma(*sptr);
        erase_maximal(sptr);
        sigma.erase(x);
        insert_simplex(sigma);
    }
}

std::size_t Toplex_map::num_simplices() const{
    return maximal_cofaces(Simplex()).size();
}

inline void Toplex_map::erase_maximal(const Simplex_ptr& sptr){
    Simplex sigma(*sptr);
    if (sptr->size()==0)
        sigma.insert(vertex_upper_bound);
    for(const Vertex& v : sigma){
        t0.at(v).erase(sptr);
        if(t0.at(v).size()==0) t0.erase(v);
    }
}

template <typename Input_vertex_range>
Vertex Toplex_map::best_index(const Input_vertex_range &vertex_range) const{
    std::size_t min = std::numeric_limits<size_t>::max();
    Vertex arg_min = vertex_upper_bound;
    for(const Vertex& v : vertex_range)
        if(!t0.count(v)) return v;
        else if(t0.at(v).size() < min)
            min = t0.at(v).size(), arg_min = v;
    return arg_min;
}

std::size_t Toplex_map::Sptr_equal::operator()(const Simplex_ptr& s1, const Simplex_ptr& s2) const {
    if (s1->size() != s2->size()) return false;
    return included(*s1,*s2);
    // inclusion tests equality for same size simplices
}

std::size_t Toplex_map::Sptr_hash::operator()(const Simplex_ptr& s) const {
    std::hash<double> h_f;
    //double hash works better than int hash
    size_t h = 0;
    for(const Vertex& v : *s)
        h += h_f(static_cast<double>(v));
    return h;
}

template <typename Input_vertex_range>
Simplex_ptr get_key(const Input_vertex_range &vertex_range){
    Simplex s(vertex_range.begin(), vertex_range.end());
    return std::make_shared<Simplex>(s);
}

template <typename Input_vertex_range1, typename Input_vertex_range2>
bool included(const Input_vertex_range1 &vertex_range1, const Input_vertex_range2 &vertex_range2){
    Simplex s2(vertex_range2.begin(), vertex_range2.end());
    for(const Vertex& v : vertex_range1)
        if(!s2.count(v)) return false;
    return true;
}

template <typename Input_vertex_range>
std::vector<Simplex> facets(const Input_vertex_range &vertex_range){
    std::vector<Simplex> facets;
    Simplex f(vertex_range.begin(), vertex_range.end());
    for(const Vertex& v : vertex_range){
        f.erase(v);
        facets.emplace_back(f);
        f.insert(v);
    }
    return facets;
}

} //namespace Gudhi

#endif /* TOPLEX_MAP_H */
