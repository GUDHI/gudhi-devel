#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include <vector>
#include <string>

#include <gudhi/Simplex_tree.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>

//typename Gudhi::Witness_complex<> Witness_complex;

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::Point_d Point_d;
typedef std::vector<Point_d> Point_Vector;
typedef CGAL::Delaunay_triangulation<K> Delaunay_triangulation;

/** \brief Write the table of the nearest landmarks to each witness
 *  to a file.
 */
template <class Value>
void write_wl( std::string file_name, std::vector< std::vector <Value> > & WL)
{
  std::ofstream ofs (file_name, std::ofstream::out);
  for (auto w : WL)
    {
      for (auto l: w)
        ofs << l << " ";
      ofs << "\n";
    }
  ofs.close();
}

/** \brief Write the coordinates of points in points to a file. 
 *
 */
void write_points( std::string file_name, std::vector< Point_d > & points)
{
  std::ofstream ofs (file_name, std::ofstream::out);
  for (auto w : points)
    {
      for (auto it = w.cartesian_begin(); it != w.cartesian_end(); ++it)
        ofs << *it << " ";
      ofs << "\n";
    }
  ofs.close();
}

/** Write edges of a witness complex in a file.
 *  The format of an edge is coordinates of u \n coordinates of v \n\n\n
 *  This format is compatible with gnuplot
 */
template< typename STree >
void write_edges(std::string file_name, STree& witness_complex, Point_Vector& landmarks)
{
  std::ofstream ofs (file_name, std::ofstream::out);
  for (auto u: witness_complex.complex_vertex_range())
    for (auto v: witness_complex.complex_vertex_range())
      {
        std::vector<int> edge = {u,v};
        if (u < v && witness_complex.find(edge) != witness_complex.null_simplex())   
          {
            for (auto it = landmarks[u].cartesian_begin(); it != landmarks[u].cartesian_end(); ++it)
              ofs << *it << " ";
            ofs << "\n";
            for (auto it = landmarks[v].cartesian_begin(); it != landmarks[v].cartesian_end(); ++it)
              ofs << *it << " ";
            ofs << "\n\n\n";
          }
      }
  ofs.close();
}

/** \brief Write triangles (tetrahedra in 3d) of a simplicial complex in a file, compatible with medit.
 *  `landmarks_ind` represents the set of landmark indices in W  
 *  `st` is the Simplex_tree to be visualized,
 *  `shr` is the Simplex_handle_range of simplices in `st` to be visualized
 *  `is2d` should be true if the simplicial complex is 2d, false if 3d
 *  `l_is_v` = landmark is vertex
 */
template <typename SimplexHandleRange,
          typename STree >
void write_witness_mesh(Point_Vector& W, std::vector<int>& landmarks_ind, STree& st, SimplexHandleRange const & shr, bool is2d, bool l_is_v, std::string file_name = "witness.mesh")
{
  std::ofstream ofs (file_name, std::ofstream::out);
  if (is2d)
    ofs << "MeshVersionFormatted 1\nDimension 2\n";
  else
    ofs << "MeshVersionFormatted 1\nDimension 3\n";
  
  if (!l_is_v)
    ofs << "Vertices\n" << W.size() << "\n";
  else
    ofs << "Vertices\n" << landmarks_ind.size() << "\n";
  
  if (l_is_v)
    for (auto p_it : landmarks_ind) {
      for (auto coord = W[p_it].cartesian_begin(); coord != W[p_it].cartesian_end() && coord != W[p_it].cartesian_begin()+3 ; ++coord)
        ofs << *coord << " ";
      ofs << "508\n";
    }
  else
    for (auto p_it : W) {
      for (auto coord = p_it.cartesian_begin(); coord != p_it.cartesian_end() && coord != p_it.cartesian_begin()+3 ; ++coord)
        ofs << *coord << " ";
      ofs << "508\n";
    }

  //  int num_triangles = W.size(), num_tetrahedra = 0;
  int num_edges = 0, num_triangles = 0, num_tetrahedra = 0;
  if (!l_is_v) {
      for (auto sh_it : shr)
        if (st.dimension(sh_it) == 1)
          num_edges++;
        else if (st.dimension(sh_it) == 2)
          num_triangles++;
        else if (st.dimension(sh_it) == 3)
          num_tetrahedra++;
      ofs << "Edges " << num_edges << "\n";
      for (auto sh_it : shr) {
           if (st.dimension(sh_it) == 1) {
            for (auto v_it : st.simplex_vertex_range(sh_it))
              ofs << landmarks_ind[v_it]+1 << " ";
           ofs << "200\n";
          }
        }
      ofs << "Triangles " << num_triangles << "\n";
      for (unsigned i = 0; i < W.size(); ++i)
        ofs << i << " " << i << " " << i << " " << "508\n";
      for (auto sh_it : shr)
        {
           if (st.dimension(sh_it) == 2) {
            for (auto v_it : st.simplex_vertex_range(sh_it))
              ofs << landmarks_ind[v_it]+1 << " ";
           ofs << "508\n";
          }
        }
      ofs << "Tetrahedra " << num_tetrahedra << "\n";
      for (auto sh_it : shr)
        {
           if (st.dimension(sh_it) == 3) {
            for (auto v_it : st.simplex_vertex_range(sh_it))
              ofs << landmarks_ind[v_it]+1 << " ";
           ofs << "250\n";
          }
        }
  }
  else {
      for (auto sh_it : shr)
        if (st.dimension(sh_it) == 1)
          num_edges++;
        else if (st.dimension(sh_it) == 2)
          num_triangles++;
        else if (st.dimension(sh_it) == 3)
          num_tetrahedra++;
      ofs << "Edges " << num_edges << "\n";
      for (auto sh_it : shr) {
           if (st.dimension(sh_it) == 1) {
            for (auto v_it : st.simplex_vertex_range(sh_it))
              ofs << v_it+1 << " ";
           ofs << "200\n";
          }
        }
      ofs << "Triangles " << num_triangles << "\n";
      for (auto sh_it : shr)
        {
           if (st.dimension(sh_it) == 2) {
            for (auto v_it : st.simplex_vertex_range(sh_it))
              ofs << v_it+1 << " ";
           ofs << "508\n";
          }
        }
      ofs << "Tetrahedra " << num_tetrahedra << "\n";
      for (auto sh_it : shr)
        {
           if (st.dimension(sh_it) == 3) {
            for (auto v_it : st.simplex_vertex_range(sh_it))
              ofs << v_it+1 << " ";
           ofs << "250\n";
          }
        }
  }
    
  ofs << "End\n";
  /*
  else
    {
      ofs << "Tetrahedra " << t.number_of_finite_full_cells()+1 << "\n";
      for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it)
        {
          if (t.is_infinite(fc_it))
            continue;
          for (auto vh_it = fc_it->vertices_begin(); vh_it != fc_it->vertices_end(); ++vh_it)
            ofs << index_of_vertex[*vh_it] << " ";
          ofs << "508\n";
        }
      ofs << nbV << " " << nbV << " " << nbV << " " << nbV << " " << 208 << "\n";
      ofs << "End\n";
    }
  */
  ofs.close();
}

void write_witness_mesh(Point_Vector& W, std::vector<int>& landmarks_ind, Gudhi::Simplex_tree<>& st, bool is2d, bool l_is_v, std::string file_name = "witness.mesh")
{
  write_witness_mesh(W, landmarks_ind, st, st.complex_simplex_range(), is2d, l_is_v, file_name);
}

/** \brief Write triangles (tetrahedra in 3d) of a Delaunay
 *  triangulation in a file, compatible with medit.
 */
void write_delaunay_mesh(Delaunay_triangulation& t, const Point_d& p, bool is2d)
{
  std::ofstream ofs ("delaunay.mesh", std::ofstream::out);
  int nbV = t.number_of_vertices()+1;
  if (is2d)
    ofs << "MeshVersionFormatted 1\nDimension 2\n";
  else
    ofs << "MeshVersionFormatted 1\nDimension 3\n";
  ofs << "Vertices\n" << nbV << "\n";
  int ind = 1; //index of a vertex
  std::map<Delaunay_triangulation::Vertex_handle, int> index_of_vertex;
  for (auto v_it = t.vertices_begin(); v_it != t.vertices_end(); ++v_it)
    {
      if (t.is_infinite(v_it))
        continue;
      // Add maximum 3 coordinates
      for (auto coord = v_it->point().cartesian_begin(); coord != v_it->point().cartesian_end() && coord != v_it->point().cartesian_begin()+3; ++coord)
        ofs << *coord << " ";
      ofs << "508\n";
      index_of_vertex[v_it] = ind++;
    }
  for (auto coord = p.cartesian_begin(); coord != p.cartesian_end(); ++coord)
    ofs << *coord << " ";
  ofs << "208\n";
  if (is2d)
    {
      ofs << "Triangles " << t.number_of_finite_full_cells()+1 << "\n";
      for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it)
        {
          if (t.is_infinite(fc_it))
            continue;
          for (auto vh_it = fc_it->vertices_begin(); vh_it != fc_it->vertices_end(); ++vh_it)
            ofs << index_of_vertex[*vh_it] << " ";
          ofs << "508\n";
        }
      ofs << nbV << " " << nbV << " " << nbV << " " << 208 << "\n";
      ofs << "End\n";
    }
  else if (p.size() == 3)
    {
      ofs << "Tetrahedra " << t.number_of_finite_full_cells()+1 << "\n";
      for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it)
        {
          if (t.is_infinite(fc_it))
            continue;
          for (auto vh_it = fc_it->vertices_begin(); vh_it != fc_it->vertices_end(); ++vh_it)
            ofs << index_of_vertex[*vh_it] << " ";
          ofs << "508\n";
        }
      ofs << nbV << " " << nbV << " " << nbV << " " << nbV << " " << 208 << "\n";
      ofs << "End\n";
    }
  else if (p.size() == 4)
    {
      ofs << "Tetrahedra " << 5*(t.number_of_finite_full_cells())+1 << "\n";
      for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it)
        {
          if (t.is_infinite(fc_it))
            continue;
          for (auto vh_it = fc_it->vertices_begin(); vh_it != fc_it->vertices_end(); ++vh_it)
            {
              for (auto vh_it2 = fc_it->vertices_begin(); vh_it2 != fc_it->vertices_end(); ++vh_it2)
                if (vh_it != vh_it2)
                  ofs << index_of_vertex[*vh_it2] << " ";
              ofs << "508\n"; 
            }
        }
      ofs << nbV << " " << nbV << " " << nbV << " " << nbV << " " << 208 << "\n";
      ofs << "End\n";
    }
  ofs.close();
}

///////////////////////////////////////////////////////////////////////
// PRINT VECTOR
///////////////////////////////////////////////////////////////////////

template <typename T>
void print_vector(std::vector<T> v)
{
  std::cout << "[";
  if (!v.empty())
    {
      std::cout << *(v.begin());
      for (auto it = v.begin()+1; it != v.end(); ++it)
        {
          std::cout << ",";
          std::cout << *it;
        }
    }
  std::cout << "]";
}


#endif
