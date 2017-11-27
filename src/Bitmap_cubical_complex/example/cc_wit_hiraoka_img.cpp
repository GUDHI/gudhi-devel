// for persistence algorithm
#include <gudhi/reader_utils.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Euclidean_witness_complex.h>
#include <gudhi/Bottleneck.h>
#include <gudhi/choose_n_farthest_points.h>

// boost libraries
// #include <boost/gil/gil_all.hpp>
// #include <boost/mpl/vector.hpp>
// #include <boost/gil/extension/io/jpeg_dynamic_io.hpp>

// standard stuff
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstddef>
#include <CGAL/Epick_d.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Delaunay_triangulation.h>
#include "../../Witness_complex/example/generators.h"
#include "output_tikz.h"

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::Point_d Point_d;
typedef CGAL::Search_traits<
  FT, Point_d,
  typename K::Cartesian_const_iterator_d,
  typename K::Construct_cartesian_const_iterator_d> Traits_base;
typedef CGAL::Search_traits_adapter<
  std::ptrdiff_t, Point_d*, Traits_base> STraits;

/** \brief Write the coordinates of points to a file. 
 *
 */
template <typename Point_d>
void write_points(std::vector< Point_d > & points, std::string file_name)
{
  std::ofstream ofs (file_name, std::ofstream::out);
  for (auto w : points)
    {
      for (auto it = w.begin(); it != w.end(); ++it)
        ofs << *it << " ";
      ofs << "\n";
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
          typename STree,
          typename Point_Vector>
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
  if (l_is_v) {
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
            ofs << st.filtration(sh_it) << "\n";
            // ofs << "200\n";
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
            // ofs << "508\n";
            ofs << st.filtration(sh_it) << "\n";
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
            // ofs << "200\n";
            ofs << st.filtration(sh_it) << "\n";
           }
        }
      ofs << "Triangles " << num_triangles << "\n";
      for (auto sh_it : shr)
        {
           if (st.dimension(sh_it) == 2) {
            for (auto v_it : st.simplex_vertex_range(sh_it))
              ofs << v_it+1 << " ";
            // ofs << "508\n";
            ofs << st.filtration(sh_it) << "\n";
          }
        }
      ofs << "Tetrahedra " << num_tetrahedra << "\n";
      for (auto sh_it : shr)
        {
           if (st.dimension(sh_it) == 3) {
            for (auto v_it : st.simplex_vertex_range(sh_it))
              ofs << v_it+1 << " ";
            // ofs << "250\n";
            ofs << st.filtration(sh_it) << "\n";
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

template <typename Point_Vector>
void write_witness_mesh(Point_Vector& W, std::vector<int>& landmarks_ind, Gudhi::Simplex_tree<>& st, bool is2d, bool l_is_v, std::string file_name = "witness.mesh")
{
  write_witness_mesh(W, landmarks_ind, st, st.complex_simplex_range(), is2d, l_is_v, file_name);
}

/** \brief Write triangles (tetrahedra in 3d) of a Delaunay
 *  triangulation in a file, compatible with medit.
 */
template <typename Delaunay_triangulation,
          typename Point_d>
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
  std::map<typename Delaunay_triangulation::Vertex_handle, int> index_of_vertex;
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


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <class Point,
          class Point_range>
void generate_witness_in_a_cube(const Point& box, unsigned nbP, Point_range& points) {
  unsigned dim = box.size();
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
  std::vector<K::Point_d> points_in_the_cube;
  generate_points_random_box(points_in_the_cube, nbP, dim+1);
  for (auto p: points_in_the_cube) {
    typedef typename Point_range::value_type Point_d;
    Point_d p_aux;
    for (auto x: p)
      p_aux.push_back(1-std::abs(x));
    typename Point_d::iterator min_it = std::min_element(p_aux.begin(), p_aux.end());
    unsigned min_i = std::distance(p_aux.begin(), min_it);
    std::vector<double> w;
    
    // for (auto x: p)
    //   std::cout << x << " ";
    // std::cout << min_i;
    // std::cout << std::endl;
    typename Point::const_reverse_iterator box_rit = box.rbegin();
    for (unsigned i = 0;
         i != p.size(); ++i)
      if (i != min_i)
        w.push_back(*(box_rit++) + (p[i]+1)/2);
        // std::cout << p[i] << " ";
    // std::cout << std::endl;
    points.push_back(Point_d(w));
  }
}

template <class Sizes,
          class Data,
          class Point_range>
void insert_witnesses(unsigned nbW, const Sizes& sizes, const Data& data, Point_range& witnesses) {
  unsigned dim = sizes.size(),
    multiplies = data.size();
  typename Data::value_type 
    max_value = *std::max_element(data.begin(), data.end()),
    sum        = 0;
  for (unsigned i = 0; i < multiplies; i++)
    sum += 1 - data[i]/max_value;
  
  std::cout << "Sum = " << sum << std::endl;
  for (unsigned i = 0; i < multiplies; i++) {
    std::vector<unsigned> inverse_box_coordinates;
    inverse_box_coordinates.reserve(dim);
    unsigned i_aux = i;
    for (unsigned j = 0; j != dim; j++) {
      inverse_box_coordinates.push_back(i_aux % sizes[dim-j-1]);
      i_aux /= sizes[dim-j-1];
    }
    generate_witness_in_a_cube(inverse_box_coordinates, std::trunc(nbW*(1-data[i]/max_value)/sum), witnesses);
  }
  std::cout << "Wrote witnesses in witnesses.dat\n";
  write_points(witnesses, "witnesses.dat");
}

template <class Sizes,
          class Data,
          class Point_range>
void insert_witnesses_naively(unsigned nbW, const Sizes& sizes, const Data& data, Point_range& witnesses) {
  typedef typename Point_range::value_type Point_d;
  
  std::vector<K::Point_d> points_in_the_cube;
  unsigned dim = sizes.size();
  generate_points_random_box(points_in_the_cube, nbW, dim);

  for (auto p: points_in_the_cube) {
    std::vector<double> w;
    int i = 0;
    for (auto x: p) {
      w.push_back((1-2*(i%2)) * (x+1)/2*sizes[i]);
      i++;
    }
    witnesses.push_back(Point_d(w));
  }
  std::cout << "Wrote witnesses in witnesses.dat\n";
  write_points(witnesses, "witnesses.dat");
}

template <class Simplicial_complex>
void fill_filtration(Simplicial_complex& sc) {
  typedef typename Simplicial_complex::Filtration_value FV;
  for (auto sh: sc.complex_simplex_range()) {
    FV f_max = 0;
    // for (auto v: sc.simplex_vertex_range(sh)) {
    //   std::cout << v << " ";
    // }
    // std::cout << " with vertex filtrations";
    for (auto v: sc.simplex_vertex_range(sh)) {
      FV f_v = sc.filtration(sc.find(std::vector<int>({v})));
      // std::cout << "";
      f_max = std::max(f_max, f_v);
    }
    sc.assign_filtration(sh, f_max);
  }
}

// template <class Boost_grayscale_image>
// void read_grayscale(const Boost_grayscale_image& src, std::vector<double>& data) {
//   for (int y = 0; y != src.height(); ++y) {
//     typename Boost_grayscale_image::x_iterator src_it = src.row_begin(y);
//     data.push_back(*src_it);
//   }
// }

template <class Delaunay_triangulation,
          class Simplicial_complex>
void convert_delaunay(const Delaunay_triangulation& del, Simplicial_complex& sc) {
  typedef typename Delaunay_triangulation::Vertex_const_handle Vertex_handle;
  typedef std::map<Vertex_handle, int> Dictionary;

  int curr_i = 0;
  Dictionary index_of_vertex;
  for (auto del_it = del.vertices_begin(); del_it != del.vertices_end(); ++del_it) {
    //index_of_vertex.emplace(std::make_pair(*del_it, curr_i));
    if (del.is_infinite(del_it))
      continue;
    index_of_vertex[del_it] = curr_i;
    curr_i++;
  }
  for (auto face_it = del.full_cells_begin(); face_it != del.full_cells_end(); ++face_it)
    if (!del.is_infinite(face_it)) {
      std::set<int> vertices;
      for (auto v_it = face_it->vertices_begin(); v_it != face_it->vertices_end(); ++v_it) {
        vertices.insert(index_of_vertex[*v_it]);
      }
      sc.insert_simplex_and_subfaces(vertices);
    }
}

static inline std::pair<double, double> compute_root_square(std::pair<double, double> input) {
  return std::make_pair(std::sqrt(input.first), std::sqrt(input.second));
}

int main(int argc, char** argv) {
  std::cout << "This program is a playground program by Siargey Kachanovich in attempt to marry the witnesses " <<
      "with the cubical complexes. It is a copy of Pawel Dlotko's example. " <<
      "This program computes persistent homology, by using bitmap_cubical_complex class, of cubical " << 
      "complexes provided in image files" <<
      "In the lines I between 2 and D+1 there are numbers of top dimensional cells in the direction I. Let " <<
      "N denote product of the numbers in the lines between 2 and D. In the lines D+2 to D+2+N there are " <<
      "filtrations of top dimensional cells. We assume that the cells are in the lexicographical order. See " <<
      "CubicalOneSphere.txt or CubicalTwoSphere.txt for example.\n" << std::endl;

  if (argc != 4) {
    std::cerr << "Wrong number of parameters. Please provide the name of a file with a Perseus style bitmap, " <<
        "the number of witnesses and the number of landmarks in the input. The program will now terminate.\n";
    return 1;
  }
  
  typedef Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
  typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
  typedef Gudhi::Simplex_tree<> Simplex_tree;
  typedef Simplex_tree::Filtration_value Filtration_value;
  typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
  typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;
  typedef Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp> Persistent_cohomology2;
  
  // Read the image
  bool dbg = false;

  //  boost::gil::gray8c_view_t src;
  std::vector<double> data;
  unsigned dimensionOfData = 2;
  
  if (dbg) {
    std::cerr << "dimensionOfData : " << dimensionOfData << std::endl;
    getchar();
  }

  // typedef boost::mpl::vector<boost::gil::gray8_image_t, boost::gil::gray16_image_t, boost::gil::rgb8_image_t, boost::gil::rgb16_image_t> my_img_types;
  // boost::gil::any_image<my_img_types> src;
  // boost::gil::jpeg_read_image(argv[1], src);
  // gil::png_read_image(argv[1], src);
  std::vector<unsigned> sizes;
  // sizes.push_back(src.height());
  // sizes.push_back(src.width());
  // read_grayscale(boost::gil::color_converted_view<boost::gil::gray8_pixel_t>(boost::gil::const_view(src)), data);
  // sizes.reserve(dimensionOfData);
  // for (size_t i = 0; i != dimensionOfData; ++i) {
  //   unsigned size_in_this_dimension;
  //   in_file >> size_in_this_dimension;
  //   sizes.push_back(size_in_this_dimension);
  //   multiplies *= size_in_this_dimension;
  //   if (dbg) {
  //     std::cerr << "size_in_this_dimension : " << size_in_this_dimension << std::endl;
  //   }
  // }
  // for (size_t i = 0; i != multiplies; ++i) {
  //   double x;
  //   in_file >> x;
  //   data.push_back(x);
  // }

  std::ifstream in_file;
  in_file.open(argv[1]);
  std::string line;
  double x;
  bool line_length_known = false;
  sizes.push_back(0);
  sizes.push_back(0);
  while (getline(in_file, line)) {
    std::istringstream iss(line);
    while (iss >> x) {
      data.push_back(x);
      if (!line_length_known)
        sizes[0]++;
    }
    line_length_known = true;
    sizes[1]++;
  }
  in_file.close();
  
  std::cout << "Width: " << sizes[0] << " Height: " << sizes[1] << std::endl;

  clock_t start, end;
  
  Bitmap_cubical_complex b(sizes, data);

  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh1(b);
  int p = 2;
  double min_persistence = 0;

  pcoh1.init_coefficients(p);  // initializes the coefficient field for homology
  pcoh1.compute_persistent_cohomology(min_persistence);

  std::string output_file_name(argv[1]);
  output_file_name += "_persistence";

  std::size_t last_in_path = output_file_name.find_last_of("/\\");

  if (last_in_path != std::string::npos) {
    output_file_name = output_file_name.substr(last_in_path+1);
  }

  std::ofstream out(output_file_name.c_str());
  pcoh1.output_diagram(out);
  out.close();

  std::cout << "Result of the regular persistence in file: " << output_file_name << "\n";
  double max_value = 0;
  for (auto x : data)
      max_value = std::max(max_value, x);
  write_persistence_diagram(output_file_name, max_value, output_file_name+".tikz.tex");

  unsigned nbW = std::atoi(argv[2]),
    nbL = std::atoi(argv[3]);
  std::vector<Point_d> witness_range, landmark_range;
  //insert_witnesses(nbW, sizes, data, witness_range);
  insert_witnesses_naively(nbW, sizes, data, witness_range);
  Gudhi::subsampling::choose_n_farthest_points(K(),
                                               witness_range,
                                               nbL,
                                               0,
                                               std::back_inserter(landmark_range));
  write_points(landmark_range, "landmarks.dat");
  
  Simplex_tree stree;
  std::vector<int> landmarks_ind;
  // Gudhi::witness_complex::Euclidean_witness_complex<K> wit(landmark_range, witness_range);
  // wit.create_complex(stree, 0, dimensionOfData);
  // write_witness_mesh(landmark_range, landmarks_ind, stree, dimensionOfData==2, false, output_file_name+".mesh");

  CGAL::Delaunay_triangulation<K> del(dimensionOfData);
  del.insert(witness_range.begin(), witness_range.end());
  convert_delaunay(del, stree);
  std::vector<Point_d> del_vertex_range;
  for (auto v_it = del.vertices_begin(); v_it != del.vertices_end(); v_it++) {
    if (del.is_infinite(*v_it))
      continue;
    del_vertex_range.push_back(v_it->point());
  }
    
  unsigned l_i = 0;
  for (auto l_it = del.vertices_begin(); l_it != del.vertices_end(); ++l_it) {
    if (del.is_infinite(*l_it))
      continue;
    unsigned pos = 0;
    for (int j = dimensionOfData-1; j != -1; j--) {
      pos *= sizes[j];
      pos += (1-2*(j%2)) * std::trunc(l_it->point().cartesian(j));
    }
    std::cout << l_i << " ";
    for (auto x: l_it->point())
      std::cout << x << " ";
    std::cout << pos << "\n";
    stree.assign_filtration(stree.find(std::vector<int>(1,l_i)), data[pos]);
    l_i++;
  }
  fill_filtration(stree);
  stree.set_dimension(dimensionOfData);
  stree.initialize_filtration();
  // std::cout << stree << std::endl;

  // persistence 2
  Persistent_cohomology2 pcoh2(stree);

  pcoh2.init_coefficients(p);  // initializes the coefficient field for homology
  pcoh2.compute_persistent_cohomology(min_persistence);

  output_file_name += "2";
  std::ofstream out2(output_file_name.c_str());
  pcoh2.output_diagram(out2);
  out2.close();
  write_persistence_diagram(output_file_name, max_value, output_file_name+".tikz.tex");

  write_witness_mesh(del_vertex_range, landmarks_ind, stree, dimensionOfData==2, false, output_file_name+".mesh");
  std::cout << "Wrote the medit file of the filtered Delaunay complex in " << output_file_name+".mesh." << std::endl;
  // write_delaunay_mesh(del, *witness_range.begin(), dimensionOfData==2);

  double max_b_distance = 0;
  for (int dim = 0; dim <= dimensionOfData; dim++) {
    std::vector< std::pair< Filtration_value , Filtration_value > >
      intervals1 = pcoh1.intervals_in_dimension(dim),
      intervals2 = pcoh2.intervals_in_dimension(dim);
    double bottleneck_distance = Gudhi::persistence_diagram::bottleneck_distance(intervals1, intervals2);
    std::cout << "In dimension " << dim << ", bottleneck distance = " << bottleneck_distance << std::endl;
    max_b_distance = std::max(max_b_distance, bottleneck_distance);
  }
  std::cout << "Bottleneck distance is " << max_b_distance << std::endl;
  
  return 0;
  
}
