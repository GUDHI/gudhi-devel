#include <iostream>
#include <random>
#include <chrono>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Toplex_map.h>
#include <gudhi/Lazy_Toplex_map.h>

using namespace Gudhi;

typedef Toplex_map::Simplex Simplex;
typedef Toplex_map::Vertex Vertex;
typedef std::pair< Simplex_tree<>::Simplex_handle, bool > typePairSimplexBool;

class ST_wrapper {

public:
  void insert_simplex(const Simplex& tau) {
    /*std::cout << "insert_simplex - " << simplexTree.num_simplices() << " - ";
    for (auto v : tau)
      std::cout << v << ", ";
    std::cout << std::endl;
    */
    simplexTree.insert_simplex_and_subfaces(tau);
  }

  bool membership(const Simplex& tau) {
    return simplexTree.find(tau) != simplexTree.null_simplex();
  }

  Vertex contraction(const Vertex x, const Vertex y) {
    // TODO (VR): edge contraction is not yet available for Simplex_tree
    return y;
  }

  std::size_t num_maximal_simplices() {
    return simplexTree.num_simplices();
  }

private:
    Simplex_tree<> simplexTree;
    void erase_max(const Simplex& sigma) {
      if(membership(sigma))
        simplexTree.remove_maximal_simplex(simplexTree.find(sigma));
    }
};


int n = 300;

int nb_insert_simplex1 = 3000;
int nb_membership1 = 4000;
int nb_contraction = 300;
int nb_insert_simplex2 = 3000;
int nb_membership2 = 400000;

Simplex random_simplex(int n, std::size_t d){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<std::size_t> dis(1, n);
    Simplex s;
    while(s.size() < d)
        s.insert(dis(gen));
    return s;
}

std::vector<Simplex> r_vector_simplices(int n, int max_d, int m){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<std::size_t> dis(1, max_d);
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
    for(int i = 1; i<=nb_contraction; i++)
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

    if (c3 > 0)
      std::cout << c1 << "\t \t" << c2 << "\t \t" << c3 << "\t \t" << K.num_maximal_simplices() << std::endl;
    else
      std::cout << c1 << "\t \t" << c2 << "\t \tN/A\t \t" << K.num_maximal_simplices() << std::endl;
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
