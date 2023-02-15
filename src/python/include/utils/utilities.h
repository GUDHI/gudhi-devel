/*    This file is part of the MMA Library - https://gitlab.inria.fr/dloiseau/multipers - which is released under MIT.
 *    See file LICENSE for full license details.
 *    Author(s):       David Loiseaux, Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include <limits>
#include <omp.h>
#include <algorithm>


namespace utils {

constexpr bool verbose = false;

using index = unsigned int;
using value_type = double;
using filtration_type = std::vector<value_type>;
using multifiltration_type = std::vector<filtration_type>;
using dimension_type = int;
using persistence_pair = std::pair<value_type, value_type>;
using boundary_type = std::vector<index>;
using boundary_matrix = std::vector<boundary_type>;
using permutation_type = std::vector<unsigned int>;
using point_type = std::vector<value_type>;
using corner_type = std::vector<value_type>;
using corners_type = std::pair<std::vector<std::vector<value_type>>, std::vector<std::vector<value_type>>>;
using bar = std::pair<corner_type, corner_type>;
using multipers_barcode = std::vector<std::vector<double>>;
const value_type inf = std::numeric_limits<value_type>::infinity();
const value_type negInf = -1 * inf;
using interval_type = std::pair<point_type, point_type>;

template<typename T>
bool is_smaller(const std::vector<T>& x, const std::vector<T>& y)
{
    for (unsigned int i = 0; i < std::min(x.size(), y.size()); i++)
        if (x[i] > y[i]) return false;
    return true;
}
template<typename T>
bool is_greater(const std::vector<T>& x, const std::vector<T>& y)
{
    for (unsigned int i = 0; i < std::min(x.size(), y.size()); i++)
        if (x[i] < y[i]) return false;
    return true;
}




struct Bar{
    Bar() : dim(-1), birth(-1), death(-1)
    {}

    Bar(dimension_type dim, int birth, int death)
        : dim(dim), birth(birth), death(death)
    {}

    dimension_type dim;
    int birth;
    int death;
};

using barcode_type = std::vector<Bar>;

struct Diagram_point{
    Diagram_point() : dim(-1), birth(-1), death(-1)
    {}

    Diagram_point(dimension_type dim,
                  value_type birth,
                  value_type death)
        : dim(dim), birth(birth), death(death)
    {}

    dimension_type dim;
    value_type birth;
    value_type death;
};
using diagram_type = std::vector<Diagram_point>;

struct MultiDiagram_point{
public:
    MultiDiagram_point() : dim(-1), birth({}), death({})
    {}

    MultiDiagram_point(dimension_type dim,
                  corner_type birth,
                  corner_type death)
    : dim(dim), birth(birth), death(death)
    {}
    dimension_type get_dimension() const {return dim;}
    corner_type get_birth() const {return birth;}
    corner_type get_death() const {return death;}
private:
    dimension_type dim;
    corner_type birth;
    corner_type death;

};
struct MultiDiagram{ // for python interface
public:
    using iterator = std::vector<MultiDiagram_point>::const_iterator;
    MultiDiagram(){}
    MultiDiagram(std::vector<MultiDiagram_point>& m) : multiDiagram(m)
    {}
    std::vector<bar> get_points(const dimension_type dimension = -1) const{ // dump for python interface
        std::vector<bar> out;
        out.reserve(multiDiagram.size());
        for (const MultiDiagram_point &pt : multiDiagram){
            if (dimension == -1 || pt.get_dimension() == dimension){
                if (pt.get_birth().size() > 0 && pt.get_death().size() > 0 && pt.get_birth()[0] != utils::inf )
					out.push_back({pt.get_birth(), pt.get_death()});
            }
        }
        out.shrink_to_fit();
        return out;
    }
    std::vector<std::vector<double>> to_multipers(const dimension_type dimension = -1) const{ // dump for python interface
        std::vector<std::vector<double>> out;
        out.reserve(multiDiagram.size());
        for (const MultiDiagram_point &pt : multiDiagram){
            if (pt.get_dimension() == dimension){
                const auto &b = pt.get_birth();
                const auto &d = pt.get_death();
                out.push_back({b[0], d[0],b[1], d[1]});
            }
        }
        out.shrink_to_fit();
        return out;
    }
    iterator begin() const {return this->multiDiagram.begin();}
    iterator end() const {return this->multiDiagram.end();}
    unsigned int size() {return this->multiDiagram.size();}
    void set(std::vector<MultiDiagram_point>& m) {this-> multiDiagram.swap(m);}
    std::vector<MultiDiagram_point>& getref() {return this->multiDiagram;}
    MultiDiagram_point& operator[](unsigned int i) {return this->multiDiagram[i];}
    MultiDiagram_point& at(const unsigned int i){
        return multiDiagram[i];
    }
private:
    std::vector<MultiDiagram_point> multiDiagram;

};
struct MultiDiagrams{
public:
    using iterator = std::vector<MultiDiagram>::const_iterator;
    using nciterator = std::vector<MultiDiagram>::iterator;
    MultiDiagrams() {}
    MultiDiagrams(unsigned int size) : multiDiagrams(size) {}
    std::vector<std::vector<std::vector<double>>> to_multipers(){
        unsigned int nsummands = this->multiDiagrams.front().size();
        unsigned int nlines = this->multiDiagrams.size();
        // std::vector<std::vector<std::vector<double>>> out(nsummands, std::vector<std::vector<double>>(nlines, std::vector<double>(5)));
        std::vector<std::vector<std::vector<double>>> out(nsummands);
        for (unsigned int i = 0; i < nsummands; i++){
            out[i].reserve(nlines);
            for(unsigned int j = 0; j < nlines; j++){
                const MultiDiagram_point &pt = this->multiDiagrams[j][i];
                if(_is_inf(pt.get_birth()) || _is_negInf(pt.get_death()))
                    out[i].push_back({0, 0, 0, 0,static_cast<double>(j)});
                else
                    out[i].push_back({pt.get_birth()[0], pt.get_death()[0], pt.get_birth()[1], pt.get_death()[1],static_cast<double>(j)});
            }
            out[i].shrink_to_fit();
        }
        return out;
    }
    using __for_python_plot_type = std::pair<std::vector<std::pair<value_type, value_type>>,std::vector<unsigned int>>;
    __for_python_plot_type _for_python_plot(dimension_type dimension=-1, double min_persistence=0){
        __for_python_plot_type out;
        auto& bars = out.first;
        auto& summand_idx= out.second;
        bars.reserve(this->multiDiagrams.size() * this->multiDiagrams[0].size()*2);
        summand_idx.reserve(this->multiDiagrams.size() * this->multiDiagrams[0].size());
        for (const MultiDiagram& multiDiagram : this->multiDiagrams){
            unsigned int count = 0;
            for (const MultiDiagram_point& bar : multiDiagram){
                const auto& birth = bar.get_birth();
                const auto& death = bar.get_death();
                if ( (dimension == -1  || bar.get_dimension() == dimension) &&  (!_is_inf(birth) && (death[0] > birth[0] + min_persistence)) ) {
                    bars.push_back({birth[0], death[0]});
                    bars.push_back({birth[1], death[1]});
                    summand_idx.push_back(count);
                }
                count++;
            }
        }
        return out;
    }
    MultiDiagram& operator[](const unsigned int i){
        return multiDiagrams[i];
    }
    MultiDiagram& at(const unsigned int i){
        return multiDiagrams[i];
    }
    iterator begin() const {return this->multiDiagrams.begin();} // cython bug : iterators like bc in bcs crash)
    iterator end() const {return this->multiDiagrams.end();}

	using barcodes = std::vector<std::vector<bar>>;
	barcodes get_points(){
		unsigned int nsummands = this->multiDiagrams.front().size();
        unsigned int nlines = this->multiDiagrams.size();
        // std::vector<std::vector<std::vector<double>>> out(nsummands, std::vector<std::vector<double>>(nlines, std::vector<double>(5)));
        barcodes out(nlines, std::vector<bar>(nsummands));
        for (unsigned int i = 0; i < nlines; i++){
            for(unsigned int j = 0; j < nsummands; j++){
                const MultiDiagram_point &pt = this->multiDiagrams[i][j];
				out[i][j] = {pt.get_birth(), pt.get_death()};
            }
        }
        return out;
	}

    unsigned int size() const {return this->multiDiagrams.size();}
private:
    std::vector<MultiDiagram> multiDiagrams;
    inline bool _is_inf(const std::vector<value_type> &truc) const{
        for (const auto coord : truc)
            if (coord != inf)   return false;
        return true;
    }
    inline bool _is_negInf(const std::vector<value_type> &truc) const{
        for (const auto coord : truc)
            if (coord != negInf)   return false;
        return true;
    }
};
} //namespace utils




namespace std{
    template<typename T>
    void sort(vector<T> &to_sort, bool(*compare)(const T&, const T&)){
        std::sort(to_sort.begin(), to_sort.end(), compare);
    }
    template<typename T>
    std::vector<T>& operator-=(std::vector<T> &result, const std::vector<T> &to_substract){
        std::transform(result.begin(), result.end(), to_substract.begin(),result.begin(), std::minus<T>());
        return result;
    }
    template<typename T>
    std::vector<T>& operator+=(std::vector<T> &result, const std::vector<T> &to_add){
        std::transform(result.begin(), result.end(), to_add.begin(),result.begin(), std::plus<T>());
        return result;
    }

} // namespace std


#endif // UTILITIES_H
