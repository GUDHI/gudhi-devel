#ifdef GUDHI_DEBUG
#define GUDHI_COX_OUTPUT_TO_HTML

#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <string>

#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vector) {
  os << "(";  
  if (vector.empty()) {
    os << ")";
    return os;
  }
  auto v_it = vector.begin();
  os << *v_it++;
  for (; v_it != vector.end(); ++v_it)
    os << ", " << *v_it;
  os << ")";
  return os;
}

/* A class to make the vector horizontal instead of vertical */
struct Straighten {
  Straighten (const Eigen::VectorXd& vector) : vector_(vector) {}
  const Eigen::VectorXd& vector_;
};

std::ostream& operator<<(std::ostream& os, const Straighten& str) {
  std::size_t size = str.vector_.size();
  os << "(" << str.vector_(0);
  if (size == 0) {
    os << ")";
    return os;
  }
  for (std::size_t i = 1; i < size; ++i)
    os << ", " << str.vector_(i);
  os << ")";
  return os;
}

struct MT_inserted_info {
  std::string qr_face_, init_face_, qr_intersection_;
  bool qr_success_, is_boundary_;
  template <class Query_result,
	    class Simplex_handle>
  MT_inserted_info(const Query_result& qr, const Simplex_handle& face, bool is_boundary)
    : qr_success_(qr.success), is_boundary_(is_boundary) {
    std::ostringstream qr_face_oss, qr_intersection_oss, init_face_oss;
    qr_face_oss << qr.face; qr_face_ = qr_face_oss.str();
    init_face_oss << face; init_face_ = init_face_oss.str();
    qr_intersection_oss << qr.intersection; qr_intersection_ = qr_intersection_oss.str();
  }
};
std::list<MT_inserted_info> mt_seed_inserted_list;

std::ostringstream mt_ostream, vis_ostream;
std::vector<std::ostringstream> cc_summary_ostream, cc_traces_ostream;

void write_head(std::ofstream& ofs) {
  ofs << "  <head>\n"
      << "    <style>\n"
      << "      span.boundary {\n"
      << "	  color: darkred;\n"
      << "	  background-color: lightgreen\n"
      << "      }\n"
      << "      span.interior {\n"
      << "	  color: navy;\n"
      << "	  background-color: yellow\n"
      << "      }\n"
      << "      ul.nav {\n"
      << "	  list-style-type: none;\n"
      << "	  margin: 0;\n"
      << "	  padding: 0;\n"
      << "	  overflow: hidden;\n"
      << "	  background-color: #333;\n"
      << "	  position: fixed;\n"
      << "	  top: 0;\n"
      << "	  width: 100%;\n"
      << "      }\n"
      << "      ul.nav li {\n"
      << "	  float: left;\n"
      << "      }\n"
      << "      ul.nav li a {\n"
      << "	  display: block;\n"
      << "	  color: white;\n"
      << "	  text-align: center;\n"
      << "	  padding: 14px 16px;\n"
      << "	  text-decoration: none;\n"
      << "      }\n"
      << "      .active {\n"
      << "	  background-color: #4CAF50;\n"
      << "      }\n"
      << "    </style>\n"
      << "  </head>\n";
}

std::string simplex_format(const std::string& simplex, bool is_boundary) {
  return (std::string)"<span class=""" + (is_boundary? "boundary": "interior") + """>" +
    (is_boundary? "B": "I") + simplex + "</span>";
}

void write_to_html(std::string file_name) {
  std::ofstream ofs(file_name + ".html", std::ofstream::out);
  ofs << "<!DOCTYPE html>\n"
      << "<html>\n";
  write_head(ofs);
  ofs << "  <body>\n"
      << "    <div style=""padding:20px;margin-top:30px;background-color:#1abc9c;"">\n"
      << "    <ul class=""nav"">\n"
      << "      <li><a href=""#mant"">Manifold tracing</a></li>\n"
      << "      <li><a href=""#cell"">Cell complex</a></li>\n"
      << "      <li><a href=""#visu"">Visualization details</a></li>\n"
      << "    </ul>\n"
      << "    </div>\n";
  ofs << "    <h1> Debug traces for " << file_name << " </h1>\n";

  ofs << "    <div id=""#mant"">\n";
  ofs << "      <h2> Manifold debug trace </h2>\n";
  ofs << "      <h3> Simplices inserted during the seed phase </h3>\n";
  ofs << "      <ul>\n";
  for (const MT_inserted_info& mt_info: mt_seed_inserted_list) {
    if (mt_info.qr_success_)
      ofs << "        <li>Inserted " << simplex_format(mt_info.qr_face_, mt_info.is_boundary_)
	  << " (initially " << simplex_format(mt_info.init_face_, mt_info.is_boundary_)
	  << ") intersection point is " << mt_info.qr_intersection_ << "</li>\n";
    else
      ofs << "        <li>Failed to insert " << mt_info.init_face_
	  << ")</li>\n";
  }
  ofs << "      </ul>\n";
  ofs << "    </div>\n";

  ofs << "    <div id=""#cell"">\n"
      << "      <h2> Cell complex debug trace </h2>\n"
      << "      <p>Go to:</p>\n"
      << "      <ul>\n"
      << "	<li><a href=""#dim0"">Dimension 0</a></li>\n"
      << "	<li><a href=""#dim1"">Dimension 1</a></li>\n"
      << "	<li><a href=""#dim2"">Dimension 2</a></li>\n"
      << "      </ul>\n";
  ofs << Straighten(Eigen::VectorXd::Random(4)) << "\n";
  ofs << "    <div id=""#visu"">\n"
      << "      <h2> Visualization details debug trace </h2>\n"
      << "    </div>\n";
  ofs << "  </body>\n";
  ofs << "</html>\n";

  ofs.close();
}

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif

