#include <iostream>
#include <random>
#include <chrono>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Lazy_Toplex_map.h>

using namespace Gudhi;

typedef Toplex_map::Vertex Vertex;
typedef Toplex_map::Simplex Simplex;

class ST_wrapper {

public:
    void insert_simplex(const Simplex& tau);
    bool membership(const Simplex& tau);
    Vertex contraction(const Vertex x, const Vertex y);
    std::size_t num_simplices();

private:
    Simplex_tree<> simplexTree;
    void erase_max(const Simplex& sigma);
};

void ST_wrapper::insert_simplex(const Simplex& tau){
    simplexTree.insert_simplex_and_subfaces(tau);
}

bool ST_wrapper::membership(const Simplex& tau) {
    return simplexTree.find(tau) != simplexTree.null_simplex();
}

void ST_wrapper::erase_max(const Simplex& sigma){
    if(membership(sigma))
        simplexTree.remove_maximal_simplex(simplexTree.find(sigma));
}

Vertex ST_wrapper::contraction(const Vertex x, const Vertex y){
    Simplex sx; sx.insert(x);
    auto hx = simplexTree.find(sx);
    if(hx != simplexTree.null_simplex())
        for(auto h : simplexTree.cofaces_simplex_range(hx,0)){
            auto sr = simplexTree.simplex_vertex_range(h);
            Simplex sigma(sr.begin(),sr.end());
            erase_max(sigma);
            sigma.erase(x);
            sigma.insert(y);
            insert_simplex(sigma);
        }
    return y;
}

std::size_t ST_wrapper::num_simplices(){
    return simplexTree.num_simplices();
}



int n = 300;

int nb_insert_simplex1 = 3000;
int nb_membership1 = 4000;
int nb_contraction = 300;
int nb_insert_simplex2 = 3000;
int nb_membership2 = 400000;

Simplex random_simplex(int n, int d){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, n);
    Simplex s;
    while(s.size()!=d)
        s.insert(dis(gen));
    return s;
}

std::vector<Simplex> r_vector_simplices(int n, int max_d, int m){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, max_d);
    std::vector<Simplex> v;
    for(int i=0; i<m; i++)
        v.push_back(random_simplex(n,dis(gen)));
    return v;
}

template<typename complex_type>
void chrono(int n, int d){
    complex_type K;
    std::vector<Simplex> simplices_insert_simplex1 = r_vector_simplices(n,d,nb_insert_simplex1);
    std::vector<Simplex> simplices_membership1 = r_vector_simplices(n,d,nb_membership1);
    std::vector<Simplex> simplices_insert_simplex2 = r_vector_simplices(n - 2*nb_contraction,d,nb_insert_simplex2);
    std::vector<Simplex> simplices_membership2 = r_vector_simplices(n - 2*nb_contraction,d,nb_membership2);
    std::chrono::time_point<std::chrono::system_clock> start, end;

    for(const Simplex& s : simplices_insert_simplex1)
        K.insert_simplex(s);

    for(const Simplex& s : simplices_membership1)
        K.membership(s);

    start = std::chrono::system_clock::now();
    for(int i = 0; i<=nb_contraction; i++)
        K.contraction(n-2*i,n-2*i-1);
    end = std::chrono::system_clock::now();
    auto c3 = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    start = std::chrono::system_clock::now();
    for(const Simplex& s : simplices_insert_simplex2)
        K.insert_simplex(s);
    end = std::chrono::system_clock::now();
    auto c1 = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    start = std::chrono::system_clock::now();
    for(const Simplex& s : simplices_membership2)
        K.membership(s);
    end = std::chrono::system_clock::now();
    auto c2 = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    std::cout << c1 << "\t \t" << c2 << "\t \t" << c3 << "\t \t" << K.num_simplices() << std::endl;
}

int main(){
    for(int d=5;d<=40;d+=5){
        std::cout << "d=" << d << " \t  Insertions \t   Membership \t Contractions \t        Size" << std::endl;
        std::cout << "T Map \t \t";
        chrono<Toplex_map>(n,d);
        std::cout << "Lazy \t \t";
        chrono<Lazy_Toplex_map>(n,d);
        if(d<=15){
            std::cout << "ST \t \t";
            chrono<ST_wrapper>(n,d);
        }
        std::cout << std::endl;
    }
}
