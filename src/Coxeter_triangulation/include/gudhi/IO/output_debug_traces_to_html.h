#ifdef GUDHI_DEBUG
#define GUDHI_COX_OUTPUT_TO_HTML

#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <regex>

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

std::string id_from_simplex(const std::string& simplex) {
  std::regex r("\\s+");
  return std::regex_replace(simplex, r, "");
}

template <typename T>
std::string to_string(const T& t) {
  std::ostringstream oss;
  oss << t;
  return oss.str();
}

struct MT_inserted_info {
  std::string qr_face_, init_face_, qr_intersection_;
  bool qr_success_, is_boundary_;
  template <class Query_result,
	    class Simplex_handle>
  MT_inserted_info(const Query_result& qr, const Simplex_handle& face, bool is_boundary)
    : qr_face_(to_string(qr.face)), init_face_(to_string(face)),
      qr_intersection_(to_string(qr.intersection)),
      qr_success_(qr.success), is_boundary_(is_boundary) {}
};
std::list<MT_inserted_info> mt_seed_inserted_list, mt_inserted_list;

struct CC_summary_info {
  std::string face_, cell_;
  template <class SC_pair>
  CC_summary_info(const SC_pair& sc_pair)
    : face_(to_string(sc_pair.first)), cell_(to_string(sc_pair.second)) {}
};
using CC_summary_list = std::list<CC_summary_info>;
std::vector<CC_summary_list> cc_interior_summary_lists, cc_boundary_summary_lists;

struct CC_detail_info {
  enum class Result_type {self, face, coface, inserted};
  std::string simplex_, trigger_;
  Result_type status_;
  std::list<std::string> faces_, post_faces_, cofaces_;
  template <class Simplex_handle>
  CC_detail_info(const Simplex_handle& simplex)
    : simplex_(to_string(simplex)) {}
};
using CC_detail_list = std::list<CC_detail_info>;
std::vector<CC_detail_list> cc_interior_detail_lists, cc_boundary_detail_lists;


std::ostringstream mt_ostream, vis_ostream;
std::vector<std::ostringstream> cc_summary_ostream, cc_traces_ostream;

void write_head(std::ofstream& ofs) {
  ofs << "  <head>\n"
      << "    <style>\n"
      << "      a.boundary {\n"
      << "	  color: darkred;\n"
      << "	  background-color: lightgreen\n"
      << "      }\n"
      << "      a.interior {\n"
      << "	  color: navy;\n"
      << "	  background-color: yellow\n"
      << "      }\n"
      << "      ul.nav {\n"
      << "	  list-style-type: none;\n"
      << "	  margin: 0;\n"
      << "	  padding: 0;\n"
      << "	  overflow: auto;\n"
      << "	  background-color: #333;\n"
      << "	  position: fixed;\n"
      // << "	  top: 0;\n"
      << "	  height: 100%;\n"
      << "	  width: 15%;\n"
      << "      }\n"
      // << "      ul.nav li {\n"
      // << "	  float: left;\n"
      // << "      }\n"
      << "      ul.nav li a {\n"
      << "	  display: block;\n"
      << "	  color: white;\n"
      << "	  text-align: left;\n"
      << "	  padding: 14px 16px;\n"
      << "	  text-decoration: none;\n"
      << "      }\n"
      << "      .active {\n"
      << "	  background-color: #4CAF50;\n"
      << "      }\n"
      << "      div {\n"
      << "        margin-left: 15%;\n"
      << "        padding: 1px 16px\n"
      << "      }\n"
      << "      div.navi {\n"
      << "        margin-left: 0%;\n"
      << "        padding: 0px 0px\n"
      << "      }\n"
      << "      h1 {\n"
      << "        margin-left: 15%;\n"
      << "        padding: 1px 16px\n"
      << "      }\n"
      << "    </style>\n"
      << "  </head>\n";
}

std::string simplex_format(const std::string& simplex, bool is_boundary) {
  std::string b_simplex = (is_boundary? "B": "I") + simplex;
  return (std::string)"<a class=""" + (is_boundary? "boundary": "interior")
    + """ href=""#" + id_from_simplex(b_simplex) + """>" + b_simplex + "</a>";
}

void write_to_html(std::string file_name) {
  std::ofstream ofs(file_name + ".html", std::ofstream::out);
  ofs << "<!DOCTYPE html>\n"
      << "<html>\n";
  write_head(ofs);
  ofs << "  <body>\n"
      << "    <div class=""navi"" style=""margin-top:30px;background-color:#1abc9c;"">\n"
      << "    <ul class=""nav"">\n"
      << "      <li><a href=""#mant"">Manifold tracing</a></li>\n"
      << "      <li><a href=""#cell"">Cell complex</a>\n"
      << "        <ul>\n";
  for (std::size_t i = 0; i < cc_interior_summary_lists.size(); ++i) {
    ofs << "          <li><a href=""#dim" << i << """>Dimension " << i << "</a>\n";
    ofs << "            <ul>\n";
    ofs << "              <li><a href=""#dim" << i << "i"">Interior</a></li>\n";
    if (i < cc_boundary_summary_lists.size()) {
      ofs << "              <li><a href=""#dim" << i << "b"">Boundary</a></li>\n";
    }
    ofs << "            </ul>\n";
    ofs << "          </li>\n";
  }
  ofs << "        </ul>\n"
      << "      </li>\n"
      << "      <li><a href=""#visu"">Visualization details</a></li>\n"
      << "    </ul>\n"
      << "    </div>\n";
  ofs << "    <h1> Debug traces for " << file_name << " </h1>\n";

  ofs << "    <div id=""mant"">\n";
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
  ofs << "      <h3> Simplices inserted during the while loop phase </h3>\n";
  ofs << "      <ul>\n";
  for (const MT_inserted_info& mt_info: mt_inserted_list) {
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

  ofs << "    <div id=""cell"">\n"
      << "      <h2> Cell complex debug trace </h2>\n"
      << "      <p>Go to:</p>\n"
      << "      <ul>\n";
  for (std::size_t i = 0; i < cc_interior_summary_lists.size(); ++i) {
    ofs << "	<li><a href=""#dim" << i << """>Dimension " << i << "</a></li>\n";
  }
  ofs << "      </ul>\n";
  for (std::size_t i = 0; i < cc_interior_summary_lists.size(); ++i) {
    ofs << "      <h3 id=""dim" << i << """> Dimension " << i << "</h3>\n";
    ofs << "      <h4 id=""dim" << i << "i""> Summary for interior simplices</h4>\n";
    if (i < cc_boundary_summary_lists.size())
      ofs << "      <p><a href=""#dim" << i << "b"">Go to boundary</a></p>\n";
    ofs << "        <ul>\n";
    for (const CC_summary_info& cc_info: cc_interior_summary_lists[i])
      ofs << "          <li id = """ << id_from_simplex("I" + cc_info.face_) << """>"
	  << simplex_format(cc_info.face_, false)
	  << " cell =" << cc_info.cell_ << "</li>\n";
    ofs << "        </ul>\n";
    ofs << "      <h4> Details for interior simplices</h4>\n";
    ofs << "        <ul>\n";
    for (const CC_detail_info& cc_info: cc_interior_detail_lists[i]) {
      ofs << "          <li> Insert_cell called for " << simplex_format(cc_info.simplex_, false)
	  << "\n";
      ofs << "            <ul>\n";
      for (const std::string& cof: cc_info.faces_)
	ofs << "              <li>Checking if " << simplex_format(cc_info.simplex_, false)
	    << " is a face of " << simplex_format(cof, false) << "\n";
      ofs << "            </ul>\n";
      ofs << "            <ul>\n";
      if (cc_info.status_ == CC_detail_info::Result_type::self) {
	ofs << "            <p><span style=""color:blue"">The simplex "
	    << simplex_format(cc_info.simplex_, false)
	    << " already exists in the cell complex!</span></p>\n";
      }
      if (cc_info.status_ == CC_detail_info::Result_type::face) {
	ofs << "            <p><span style=""color:red"">The simplex "
	    << simplex_format(cc_info.simplex_, false) << " is a face of the simplex "
	    << simplex_format(cc_info.trigger_, false) << "!</span><br>\n";
	ofs << "              <ul>\n";
	for (const std::string post_face: cc_info.post_faces_)
	  ofs << "                <li id = """ << id_from_simplex("I" + post_face) << """>"
	      << "Post deleting " << simplex_format(post_face, false) << "</li>\n";
	ofs << "              </ul>\n";
	ofs << "            </p>\n";
	ofs << "            <p id = """ << id_from_simplex("I" + cc_info.trigger_) << """>"
	    << "Deleting " << simplex_format(cc_info.trigger_, false) << "</p>\n";
      }
      for (const std::string& fac: cc_info.cofaces_)
	ofs << "              <li>Checking if " << simplex_format(cc_info.simplex_, false)
	    << " is a coface of " << simplex_format(fac, false) << "\n";
      if (cc_info.status_ == CC_detail_info::Result_type::coface) {
	ofs << "            <p><span style=""color:darkorange"">The simplex "
	    << simplex_format(cc_info.simplex_, false) << " is a coface of the simplex "
	    << simplex_format(cc_info.trigger_, false) << "!</span><p>\n";
      }
      if (cc_info.status_ == CC_detail_info::Result_type::inserted) {
	ofs << "            <p><span style=""color:green"">Successfully inserted "
	    << simplex_format(cc_info.simplex_, false) << "!</span><p>\n";
      }
      ofs << "            </ul>\n";
      ofs << "          </li>\n";
    }
    ofs << "        </ul>\n";

    if (i < cc_boundary_summary_lists.size()) {
      ofs << "      <h4 id=""dim" << i << "b""> Summary for boundary simplices</h4>\n";
      ofs << "      <p><a href=""#dim" << i << "i"">Go to interior</a></p>\n";
      ofs << "        <ul>\n";
      for (const CC_summary_info& cc_info: cc_boundary_summary_lists[i])
	ofs << "          <li  id = """ << id_from_simplex("B" + cc_info.face_) << """>"
	    << simplex_format(cc_info.face_, true)
	    << " cell =" << cc_info.cell_ << "</li>\n";
      ofs << "        </ul>\n";
    }
    if (i < cc_boundary_detail_lists.size()) {
      ofs << "      <h4> Details for boundary simplices</h4>\n"
	  << "        <ul>\n";
      for (const CC_detail_info& cc_info: cc_boundary_detail_lists[i]) {
	ofs << "          <li>" << simplex_format(cc_info.simplex_, true);
	ofs << "            <ul>\n";
	for (const std::string& cof: cc_info.faces_)
	  ofs << "              <li>Checking if " << simplex_format(cc_info.simplex_, true)
	      << " is a face of " << simplex_format(cof, true) << "\n";
	ofs << "            </ul>\n";
	ofs << "            <ul>\n";
	if (cc_info.status_ == CC_detail_info::Result_type::self) {
	  ofs << "            <p><span style=""color:blue"">The simplex "
	      << simplex_format(cc_info.simplex_, true)
	      << " already exists in the cell complex!</span></p>\n";
	}
	if (cc_info.status_ == CC_detail_info::Result_type::face) {
	  ofs << "            <p><span style=""color:red"">The simplex "
	      << simplex_format(cc_info.simplex_, true) << " is a face of the simplex "
	      << simplex_format(cc_info.trigger_, true) << "!</span><br>\n";
	  ofs << "              <ul>\n";
	  for (const std::string post_face: cc_info.post_faces_)
	    ofs << "                <li id=""" << id_from_simplex("B" + post_face)
		<< """>Post deleting " << simplex_format(post_face, true)
		<< "</li>\n";
	  ofs << "              </ul>\n";
	  ofs << "            </p>\n";
	  ofs << "            <p id=""" << id_from_simplex(cc_info.trigger_)
	      << """>Deleting " << simplex_format(cc_info.trigger_, true) << "</p>\n";
	}
	for (const std::string& fac: cc_info.cofaces_)
	  ofs << "              <li>Checking if " << simplex_format(cc_info.simplex_, true)
	      << " is a coface of " << simplex_format(fac, true) << "\n";
	ofs << "            </ul>\n";
	ofs << "          </li>\n";
	if (cc_info.status_ == CC_detail_info::Result_type::coface) {
	  ofs << "            <p><span style=""color:darkorange"">The simplex "
	      << simplex_format(cc_info.simplex_, true) << " is a coface of the simplex "
	      << simplex_format(cc_info.trigger_, true) << "!</span><p>\n";
	}
	if (cc_info.status_ == CC_detail_info::Result_type::inserted) {
	  ofs << "            <p><span style=""color:green"">Successfully inserted "
	      << simplex_format(cc_info.simplex_, true) << "!</span><p>\n";
	}
      }
      ofs << "        </ul>\n";
    }
  }  
  ofs << "    </div>\n";
  
  ofs << "    <div id=""visu"">\n"
      << "      <h2> Visualization details debug trace </h2>\n"
      << "    </div>\n";
  ofs << "  </body>\n";
  ofs << "</html>\n";

  ofs.close();
}

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif

