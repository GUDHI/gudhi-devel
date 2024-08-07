/* Based on Ripser commit 140670f2c76997404601e43d8054151f46be9fd7
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 *      - 2024 Marc Glisse: Heavy refactoring
*/

/*

 Ripser: a lean C++ code for computation of Vietoris-Rips persistence barcodes

 MIT License

 Copyright (c) 2015–2021 Ulrich Bauer

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

 You are under no obligation whatsoever to provide any bug fixes, patches, or
 upgrades to the features, functionality or performance of the source code
 ("Enhancements") to anyone; however, if you choose to make your Enhancements
 available either publicly, or directly to the author of this software, without
 imposing a separate written license agreement for such Enhancements, then you
 hereby grant the following license: a non-exclusive, royalty-free perpetual
 license to install, use, modify, prepare derivative works, incorporate into
 other computer software, distribute, and sublicense such enhancements or
 derivative works thereof, in binary and source code form.

*/

//#define GUDHI_INDICATE_PROGRESS

#include <gudhi/ripser.h>
#include <fstream>
#include <sstream>
#include <iostream>

using Gudhi::ripser::ripser;

struct Params {
  typedef float value_t;
  typedef int8_t dimension_t;
  typedef int vertex_t;
  typedef uint_least32_t coefficient_t;
};

typedef Params::value_t value_t;
typedef Params::dimension_t dimension_t;
typedef Params::vertex_t vertex_t;
typedef Params::coefficient_t coefficient_t;

typedef Gudhi::ripser::Full_distance_matrix<Params> Full_distance_matrix;
typedef Gudhi::ripser::Compressed_distance_matrix<Params, Gudhi::ripser::LOWER_TRIANGULAR> Compressed_lower_distance_matrix;
typedef Gudhi::ripser::Compressed_distance_matrix<Params, Gudhi::ripser::UPPER_TRIANGULAR> Compressed_upper_distance_matrix;
typedef Gudhi::ripser::Sparse_distance_matrix<Params> Sparse_distance_matrix;
typedef Gudhi::ripser::Euclidean_distance_matrix<Params> Euclidean_distance_matrix;

enum file_format {
  LOWER_DISTANCE_MATRIX,
  UPPER_DISTANCE_MATRIX,
  DISTANCE_MATRIX,
  POINT_CLOUD,
  DIPHA,
  SPARSE,
  BINARY
};

static constexpr uint16_t endian_check=0xff00;
static bool is_big_endian() { return *reinterpret_cast<const uint8_t*>(&endian_check); }

template <typename T> T read(std::istream& input_stream) {
  T result;
  char* p = reinterpret_cast<char*>(&result);
  if (input_stream.read(p, sizeof(T)).gcount() != sizeof(T)) return T();
  if (is_big_endian()) std::reverse(p, p + sizeof(T));
  return result;
}

Euclidean_distance_matrix read_point_cloud(std::istream& input_stream) {
  std::vector<std::vector<value_t>> points;

  std::string line;
  value_t value;
  while (std::getline(input_stream, line)) {
    std::vector<value_t> point;
    std::istringstream s(line);
    while (s >> value) {
      point.push_back(value);
      s.ignore();
    }
    if (!point.empty()) points.push_back(point);
    assert(point.size() == points.front().size());
  }

  std::size_t d = points.front().size();
  Euclidean_distance_matrix eucl_dist(std::move(points));
  vertex_t n = eucl_dist.size();
  std::cout << "point cloud with " << n << " points in dimension "
    << d << std::endl;

  return eucl_dist;
}

Sparse_distance_matrix read_sparse_distance_matrix(std::istream& input_stream) {
  typedef typename Sparse_distance_matrix::vertex_diameter_t vertex_diameter_t;
  std::vector<std::vector<vertex_diameter_t>> neighbors;
  std::size_t num_edges = 0;

  std::string line;
  while (std::getline(input_stream, line)) {
    std::istringstream s(line);
    std::size_t i, j;
    value_t value;
    s >> i;
    s.ignore();
    s >> j;
    s.ignore();
    s >> value;
    s.ignore();
    if (i != j) {
      neighbors.resize(std::max({neighbors.size(), i + 1, j + 1}));
      neighbors[i].emplace_back(j, value);
      neighbors[j].emplace_back(i, value);
      ++num_edges;
    }
  }

  for (std::size_t i = 0; i < neighbors.size(); ++i)
    std::sort(neighbors[i].begin(), neighbors[i].end());

  return Sparse_distance_matrix(std::move(neighbors), num_edges);
}

Compressed_lower_distance_matrix read_lower_distance_matrix(std::istream& input_stream) {
  std::vector<value_t> distances;
  value_t value;
  while (input_stream >> value) {
    distances.push_back(value);
    input_stream.ignore();
  }

  return Compressed_lower_distance_matrix(std::move(distances));
}

Compressed_lower_distance_matrix read_upper_distance_matrix(std::istream& input_stream) {
  std::vector<value_t> distances;
  value_t value;
  while (input_stream >> value) {
    distances.push_back(value);
    input_stream.ignore();
  }

  return Compressed_lower_distance_matrix(Compressed_upper_distance_matrix(std::move(distances)));
}

Compressed_lower_distance_matrix read_distance_matrix(std::istream& input_stream) {
  std::vector<value_t> distances;

  std::string line;
  value_t value;
  for (int i = 0; std::getline(input_stream, line); ++i) {
    std::istringstream s(line);
    for (int j = 0; j < i && s >> value; ++j) {
      distances.push_back(value);
      s.ignore();
    }
  }

  return Compressed_lower_distance_matrix(std::move(distances));
}

Compressed_lower_distance_matrix read_dipha(std::istream& input_stream) {
  if (read<int64_t>(input_stream) != 8067171840) {
    std::cerr << "input is not a Dipha file (magic number: 8067171840)" << std::endl;
    exit(-1);
  }

  if (read<int64_t>(input_stream) != 7) {
    std::cerr << "input is not a Dipha distance matrix (file type: 7)" << std::endl;
    exit(-1);
  }

  vertex_t n = read<int64_t>(input_stream);

  std::vector<value_t> distances;

  for (vertex_t i = 0; i < n; ++i)
    for (vertex_t j = 0; j < n; ++j)
      if (i > j)
        distances.push_back(read<double>(input_stream));
      else
        read<double>(input_stream);

  return Compressed_lower_distance_matrix(std::move(distances));
}

Compressed_lower_distance_matrix read_binary(std::istream& input_stream) {
  std::vector<value_t> distances;
  while (!input_stream.eof()) distances.push_back(read<value_t>(input_stream));
  return Compressed_lower_distance_matrix(std::move(distances));
}

Compressed_lower_distance_matrix read_file(std::istream& input_stream, const file_format format) {
  switch (format) {
    case LOWER_DISTANCE_MATRIX:
      return read_lower_distance_matrix(input_stream);
    case UPPER_DISTANCE_MATRIX:
      return read_upper_distance_matrix(input_stream);
    case DISTANCE_MATRIX:
      return read_distance_matrix(input_stream);
    case POINT_CLOUD:
      return read_point_cloud(input_stream);
    case DIPHA:
      return read_dipha(input_stream);
    default:
      return read_binary(input_stream);
  }
}

void print_usage_and_exit(int exit_code) {
  std::cerr
    << "Usage: "
    << "ripser "
    << "[options] [filename]" << std::endl
    << std::endl
    << "Options:" << std::endl
    << std::endl
    << "  --help           print this screen" << std::endl
    << "  --format         use the specified file format for the input. Options are:"
    << std::endl
    << "                     lower-distance (lower triangular distance matrix)"
    << std::endl
    << "                     upper-distance (upper triangular distance matrix)" << std::endl
    << "         (default:)  distance       (distance matrix; only lower triangular part is read)" << std::endl
    << "                     point-cloud    (point cloud in Euclidean space)" << std::endl
    << "                     dipha          (distance matrix in DIPHA file format)" << std::endl
    << "                     sparse         (sparse distance matrix in sparse triplet format)"
    << std::endl
    << "                     binary         (lower triangular distance matrix in binary format)"
    << std::endl
    << "  --dim <k>        compute persistent homology up to dimension k" << std::endl
    << "  --threshold <t>  compute Rips complexes up to diameter t" << std::endl
    << "  --modulus <p>    compute homology with coefficients in the prime field Z/pZ"
    << std::endl
    << "  --ratio <r>      only show persistence pairs with death/birth ratio > r" << std::endl
    << std::endl;
  exit(exit_code);
}

int main(int argc, char** argv) {
  const char* filename = nullptr;

  file_format format = DISTANCE_MATRIX;

  dimension_t dim_max = 1;
  value_t threshold = std::numeric_limits<value_t>::max();
  float ratio = 1;
  coefficient_t modulus = 2;

  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg == "--help") {
      print_usage_and_exit(0);
    } else if (arg == "--dim") {
      std::string parameter = std::string(argv[++i]);
      std::size_t next_pos;
      dim_max = std::stol(parameter, &next_pos);
      if (next_pos != parameter.size()) print_usage_and_exit(-1);
    } else if (arg == "--threshold") {
      std::string parameter = std::string(argv[++i]);
      std::size_t next_pos;
      threshold = std::stof(parameter, &next_pos);
      if (next_pos != parameter.size()) print_usage_and_exit(-1);
    } else if (arg == "--ratio") {
      std::string parameter = std::string(argv[++i]);
      std::size_t next_pos;
      ratio = std::stof(parameter, &next_pos);
      if (next_pos != parameter.size()) print_usage_and_exit(-1);
    } else if (arg == "--format") {
      std::string parameter = std::string(argv[++i]);
      if (parameter.rfind("lower", 0) == 0)
        format = LOWER_DISTANCE_MATRIX;
      else if (parameter.rfind("upper", 0) == 0)
        format = UPPER_DISTANCE_MATRIX;
      else if (parameter.rfind("dist", 0) == 0)
        format = DISTANCE_MATRIX;
      else if (parameter.rfind("point", 0) == 0)
        format = POINT_CLOUD;
      else if (parameter == "dipha")
        format = DIPHA;
      else if (parameter == "sparse")
        format = SPARSE;
      else if (parameter == "binary")
        format = BINARY;
      else
        print_usage_and_exit(-1);
    } else if (arg == "--modulus") {
      std::string parameter = std::string(argv[++i]);
      std::size_t next_pos;
      modulus = std::stol(parameter, &next_pos);
      if (next_pos != parameter.size() || !Gudhi::ripser::is_prime(modulus)) print_usage_and_exit(-1);
    } else {
      if (filename) { print_usage_and_exit(-1); }
      filename = argv[i];
    }
  }

  std::ifstream file_stream(filename);
  if (filename && file_stream.fail()) {
    std::cerr << "couldn't open file " << filename << std::endl;
    exit(-1);
  }

  auto output_dim = [](dimension_t dim) {
    std::cout << "persistence intervals in dim " << (int)dim << ":" << std::endl;
  };
  auto output_pair = [ratio](value_t birth, value_t death) {
#ifdef GUDHI_INDICATE_PROGRESS
    // Not necessary if we redirect stdout
    std::cerr << Gudhi::ripser::clear_line << std::flush;
#endif
    if (death == std::numeric_limits<value_t>::infinity())
      std::cout << " [" << birth << ", )" << std::endl;
    else if (death > birth * ratio)
      std::cout << " [" << birth << "," << death << ")" << std::endl;
  };
  if (format == SPARSE) {
    Sparse_distance_matrix dist =
      read_sparse_distance_matrix(filename ? file_stream : std::cin);
    std::cout << "sparse distance matrix with " << dist.size() << " points and "
      << dist.num_edges << "/" << (dist.size() * (dist.size() - 1)) / 2 << " entries"
      << std::endl;

    ripser(std::move(dist), dim_max, threshold, modulus, output_dim, output_pair);
  } else if (format == POINT_CLOUD && threshold < std::numeric_limits<value_t>::max()) {
    Sparse_distance_matrix dist(read_point_cloud(filename ? file_stream : std::cin), threshold);
    ripser(std::move(dist), dim_max, threshold, modulus, output_dim, output_pair);
  } else {
    Compressed_lower_distance_matrix dist =
      read_file(filename ? file_stream : std::cin, format);

    value_t min = std::numeric_limits<value_t>::infinity(),
            max = -std::numeric_limits<value_t>::infinity(), max_finite = max;
    std::size_t num_edges = 0;

    value_t enclosing_radius = std::numeric_limits<value_t>::infinity();
    if (threshold >= std::numeric_limits<value_t>::max()) {
      for (vertex_t i = 0; i < dist.size(); ++i) {
        value_t r_i = -std::numeric_limits<value_t>::infinity();
        for (vertex_t j = 0; j < dist.size(); ++j) r_i = std::max(r_i, dist(i, j));
        enclosing_radius = std::min(enclosing_radius, r_i);
      }
    }

    for (auto d : dist.distances) {
      min = std::min(min, d);
      max = std::max(max, d);
      if (d != std::numeric_limits<value_t>::infinity()) max_finite = std::max(max_finite, d);
      if (d <= threshold) ++num_edges;
    }
    std::cout << "value range: [" << min << "," << max_finite << "]" << std::endl;

    if (threshold >= std::numeric_limits<value_t>::max()) {
      std::cout << "distance matrix with " << dist.size()
        << " points, using threshold at enclosing radius " << enclosing_radius
        << std::endl;
      ripser(std::move(dist), dim_max, enclosing_radius, modulus, output_dim, output_pair);
    } else {
      std::cout << "sparse distance matrix with " << dist.size() << " points and "
        << num_edges << "/" << (dist.size() * (dist.size() - 1)) / 2 << " entries"
        << std::endl;

      ripser(Sparse_distance_matrix(std::move(dist), threshold), dim_max, threshold, modulus, output_dim, output_pair);
    }
  }
  return 0;
}
