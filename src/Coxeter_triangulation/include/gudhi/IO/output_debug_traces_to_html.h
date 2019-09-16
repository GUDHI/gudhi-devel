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
  std::regex r("\\s+"), r2("\\(|\\)|\\{|\\}"), r3(","), r4("\\["), r5("\\]");
  std::string output = std::regex_replace(simplex, r, "");
  output = std::regex_replace(output, r2, ":");
  output = std::regex_replace(output, r3, ".");
  output = std::regex_replace(output, r4, "_");
  output = std::regex_replace(output, r5, "");
  return output;
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
  enum class Result_type {self, face, coface, inserted, join_single, join_is_face};
  std::string simplex_, trigger_, init_simplex_;
  Result_type status_;
  bool join_trigger_ = false;
  std::list<std::string> faces_, post_faces_, cofaces_;
  template <class Simplex_handle>
  CC_detail_info(const Simplex_handle& simplex)
    : simplex_(to_string(simplex)) {}
};
using CC_detail_list = std::list<CC_detail_info>;
std::vector<CC_detail_list> cc_interior_detail_lists, cc_boundary_detail_lists;
std::vector<CC_detail_list> cc_interior_insert_detail_lists, cc_boundary_insert_detail_lists;

struct CC_prejoin_info {
  enum class Result_type {join_single, join_is_face, join_different, join_same};
  std::string simplex_, join_;
  std::vector<std::string> faces_;
  std::size_t dimension_;
  Result_type status_;
  template <class Simplex_handle>
  CC_prejoin_info(const Simplex_handle& simplex)
    : simplex_(to_string(simplex)), dimension_(simplex.dimension()) {}
};
using CC_prejoin_list = std::list<CC_prejoin_info>;
std::vector<CC_prejoin_list> cc_interior_prejoin_lists, cc_boundary_prejoin_lists;


struct CC_join_info {
  enum class Result_type {self, face, coface, inserted, join_single, join_is_face};
  std::string simplex_, join_, trigger_;
  Result_type status_;
  std::list<std::string> boundary_faces_;
  std::list<std::string> faces_, post_faces_, cofaces_;
  template <class Simplex_handle>
  CC_join_info(const Simplex_handle& simplex)
    : simplex_(to_string(simplex)) {}
};
bool join_switch = false;
std::vector<CC_detail_list> cc_interior_join_detail_lists, cc_boundary_join_detail_lists;

std::map<std::string, std::string> cell_vlist_map;
std::map<std::string, std::string> simplex_vlist_map;

std::ostringstream mt_ostream, vis_ostream;
std::vector<std::ostringstream> cc_summary_ostream, cc_traces_ostream;

std::string simplex_format(const std::string& simplex, bool is_boundary) {
  std::string b_simplex = (is_boundary? "B": "I") + simplex;
  std::string tooltiptext;
  auto it = simplex_vlist_map.find(b_simplex);
  if (it == simplex_vlist_map.end())
    tooltiptext = "deleted";
  else
    tooltiptext = simplex_vlist_map.at(b_simplex);
  return (std::string)"<a class=\"" + (is_boundary? "boundary": "interior")
    + "\" href=\"#" + id_from_simplex(b_simplex) + "\">" + b_simplex
    + "<span class=\"tooltiptext\">" + tooltiptext + "</span></a>";
}

std::string simplex_format(const std::string& b_simplex) {
  bool is_boundary = b_simplex[0] == 'B';
  std::string tooltiptext;
  auto it = simplex_vlist_map.find(b_simplex);
  if (it == simplex_vlist_map.end())
    tooltiptext = "deleted";
  else
    tooltiptext = simplex_vlist_map.at(b_simplex);
  return (std::string)"<a class=\"" + (is_boundary? "boundary": "interior")
    + "\" href=\"#" + id_from_simplex(b_simplex) + "\">" + b_simplex
    + "<span class=\"tooltiptext\">" + tooltiptext + "</span></a>";
}


void write_head(std::ofstream& ofs) {
  ofs << "  <head>\n"
      << "    <title>Cell complex debug trace</title>\n"
      << "    <style>\n"
      << "      a.boundary {\n"
      << "        position: relative;\n"
      << "        display: inline-block;\n"
      << "	  color: darkred;\n"
      << "	  background-color: lightgreen\n"
      << "      }\n"
      << "      a.interior {\n"
      << "        position: relative;\n"
      << "        display: inline-block;\n"
      << "	  color: navy;\n"
      << "	  background-color: yellow\n"
      << "      }\n"
      << "      .tooltiptext {\n"
      << "        visibility: hidden;\n"
      << "        width: 120px;\n"
      << "        background-color: #555;\n"
      << "        color: #fff;\n"
      << "        text-align: center;\n"
      << "        padding: 5px 0;\n"
      << "        border-radius: 6px;\n"
      << "        position: absolute;\n"
      << "        z-index: 1;\n"
      << "        bottom: 125%;\n"
      << "        left: 50%;\n"
      << "        margin-left: -60px;\n"
      << "        opacity: 0;\n"
      << "        transition: opacity 0.3s;\n"
      << "      }\n"
      << "      .boundary .tooltiptext::after {\n"
      << "        content: \"\";\n"
      << "        position: absolute;\n"
      << "        top: 100%;\n"
      << "        left: 50%;\n"
      << "        margin-left: -5px;\n"
      << "        border-width: 5px;\n"
      << "        border-style: solid;\n"
      << "        border-color: #555 transparent transparent transparent;\n"
      << "      }\n"
      << "      .interior .tooltiptext::after {\n"
      << "        content: \"\";\n"
      << "        position: absolute;\n"
      << "        top: 100%;\n"
      << "        left: 50%;\n"
      << "        margin-left: -5px;\n"
      << "        border-width: 5px;\n"
      << "        border-style: solid;\n"
      << "        border-color: #555 transparent transparent transparent;\n"
      << "      }\n"
      << "      .boundary:hover .tooltiptext {\n"
      << "        visibility: visible;\n"
      << "        opacity: 1;\n"
      << "      }\n"    
      << "      .interior:hover .tooltiptext {\n"
      << "        visibility: visible;\n"
      << "        opacity: 1;\n"
      << "      }\n"    
      << "      ul.nav {\n"
      << "	  list-style-type: none;\n"
      << "	  margin: 0;\n"
      << "	  padding: 0;\n"
      << "	  overflow: auto;\n"
      << "	  background-color: #333;\n"
      << "	  position: fixed;\n"
      << "	  height: 100%;\n"
      << "	  width: 15%;\n"
      << "      }\n"
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

void write_nav(std::ofstream& ofs) {
  ofs << "    <div class=\"navi\" style=\"margin-top:30px;background-color:#1abc9c;\">\n"
      << "    <ul class=\"nav\">\n"
      << "      <li><a href=\"#mant\">Manifold tracing</a></li>\n"
      << "      <li><a href=\"#cell\">Cell complex</a>\n"
      << "        <ul>\n";
  for (std::size_t i = 0; i < cc_interior_summary_lists.size(); ++i) {
    ofs << "          <li><a href=\"#dim" << i << "\">Dimension " << i << "</a>\n";
    ofs << "            <ul>\n";
    ofs << "              <li><a href=\"#dim" << i << "i\">Interior</a></li>\n";
    if (i < cc_boundary_summary_lists.size()) {
      ofs << "              <li><a href=\"#dim" << i << "b\">Boundary</a></li>\n";
    }
    ofs << "            </ul>\n";
    ofs << "          </li>\n";
  }
  ofs << "        </ul>\n"
      << "      </li>\n"
      << "      <li><a href=\"#visu\">Visualization details</a></li>\n"
      << "    </ul>\n"
      << "    </div>\n";
}

void write_mt(std::ofstream& ofs) {
  ofs << "    <div id=\"mant\">\n";
  ofs << "      <h2> Manifold debug trace </h2>\n";
  ofs << "      <h3> Simplices inserted during the seed phase </h3>\n";
  ofs << "      <ul>\n";
  for (const MT_inserted_info& mt_info: mt_seed_inserted_list) {
    if (mt_info.qr_success_) {
      ofs << "        <li>Inserted " << simplex_format(mt_info.qr_face_, mt_info.is_boundary_);
      if (mt_info.qr_face_ != mt_info.init_face_)
	ofs << " (initially " << simplex_format(mt_info.init_face_, mt_info.is_boundary_) << ")";
      ofs << " intersection point is " << mt_info.qr_intersection_ << "</li>\n";
    }
    else
      ofs << "        <li>Failed to insert " << mt_info.init_face_
	  << "</li>\n";
  }
  ofs << "      </ul>\n";
  ofs << "      <h3> Simplices inserted during the while loop phase </h3>\n";
  ofs << "      <ul>\n";
  for (const MT_inserted_info& mt_info: mt_inserted_list) {
    if (mt_info.qr_success_) {
      ofs << "        <li>Inserted " << simplex_format(mt_info.qr_face_, mt_info.is_boundary_);
      if (mt_info.qr_face_ != mt_info.init_face_)
	ofs << " (initially " << simplex_format(mt_info.init_face_, mt_info.is_boundary_) << ")";
      ofs << " intersection point is " << mt_info.qr_intersection_ << "</li>\n";
    }
    else
      ofs << "        <li>Failed to insert " << mt_info.init_face_
	  << ")</li>\n";
  }
  ofs << "      </ul>\n";
  ofs << "    </div>\n";
}

void write_cc(std::ofstream& ofs) {
  ofs << "    <div id=\"cell\">\n"
      << "      <h2> Cell complex debug trace </h2>\n"
      << "      <p>Go to:</p>\n"
      << "      <ul>\n";
  for (std::size_t i = 0; i < cc_interior_summary_lists.size(); ++i) {
    ofs << "	<li><a href=\"#dim" << i << "\">Dimension " << i << "</a></li>\n";
  }
  ofs << "      </ul>\n";
  for (std::size_t i = 0; i < cc_interior_summary_lists.size(); ++i) {
    ofs << "      <h3 id=\"dim" << i << "\"> Dimension " << i << "</h3>\n";
    ofs << "      <h4 id=\"dim" << i << "i\"> Summary for interior simplices</h4>\n";
    if (i < cc_boundary_summary_lists.size())
      ofs << "      <p><a href=\"#dim" << i << "b\">Go to boundary</a></p>\n";
    ofs << "        <ul>\n";
    for (const CC_summary_info& cc_info: cc_interior_summary_lists[i])
      ofs << "          <li id = \"" << id_from_simplex("I" + cc_info.face_) << "\">"
	  << simplex_format(cc_info.face_, false)
	  << " cell =" << cc_info.cell_ << "</li>\n";
    ofs << "        </ul>\n";
    ofs << "      <h4> Prejoin state of the interior cells of dimension " << i << "</h4>\n";
    auto prejoin_it = cc_interior_prejoin_lists[i].begin();
    while (prejoin_it != cc_interior_prejoin_lists[i].end()) {
      std::size_t j = prejoin_it->dimension_;
      ofs << "        <h5>" << j << "-dimensional ambient simplices</h5>\n";
      ofs << "          <ul>\n";
      for (; prejoin_it->dimension_ == j; ++prejoin_it) {
	ofs << "          <li>" << simplex_format(prejoin_it->simplex_, false)
	    << " join = " << simplex_format(prejoin_it->join_, false)
	    << " boundary:\n"
	    << "            <ul>\n";
	for (const auto& face: prejoin_it->faces_)
	  ofs << "              <li>" << simplex_format(face) << "</li>";
	ofs << "            </ul>\n";
	switch (prejoin_it->status_) {
	case (CC_prejoin_info::Result_type::join_single):
	  ofs << "            <p style=\"color: red\">Deleted "
	      << simplex_format(prejoin_it->simplex_, false)
	      << " as it has a single face.</p>";
	  break;
	case (CC_prejoin_info::Result_type::join_is_face):
	  ofs << "            <p style=\"color: red\">Deleted "
	      << simplex_format(prejoin_it->simplex_, false)
	      << " as its join " << simplex_format(prejoin_it->join_, false)
	      << " is one of the faces.</p>";
	  break;
	case (CC_prejoin_info::Result_type::join_different):
	  ofs << "            <p style=\"color: magenta\">Deleted " << simplex_format(prejoin_it->simplex_, false)
	      << " and replaced by its join " << simplex_format(prejoin_it->join_, false)
	      << ".</p>";
	  break;
	case (CC_prejoin_info::Result_type::join_same):
	  ofs << "            <p style=\"color: green\">Kept " << simplex_format(prejoin_it->simplex_, false)
	      << ".</p>";
	}
	ofs << "          </li>";
      }
      ofs << "          </ul>\n";
    }
    ofs << "      <h4> Details for interior simplices</h4>\n";
    ofs << "        <ul>\n";
    for (const CC_detail_info& cc_info: cc_interior_detail_lists[i]) {
      if (cc_info.status_ == CC_detail_info::Result_type::join_single) {
	ofs << "          <li style=\"color:magenta\" id = \""
	    << id_from_simplex("I" + cc_info.simplex_) << "\"> Simplex "
	    << simplex_format(cc_info.simplex_, false) << " has only one face ("
	    << simplex_format(cc_info.trigger_, false) << ") and is deleted.";
	continue;
      }
      if (cc_info.status_ == CC_detail_info::Result_type::join_single) {
	ofs << "          <li style=\"color:darkmagenta\" id = \""
	    << id_from_simplex("I" + cc_info.simplex_) << "\"> The join of the simplex "
	    << simplex_format(cc_info.simplex_, false) << " is one of its faces ("
	    << simplex_format(cc_info.trigger_, false) << "), hence it is is deleted.";
	continue;
      }
      ofs << "          <li> Insert_cell called for " << simplex_format(cc_info.simplex_, false)
	  << "\n";
      ofs << "            <ul>\n";
      // for (const std::string& cof: cc_info.faces_)
      // 	ofs << "              <li>Checking if " << simplex_format(cc_info.simplex_, false)
      // 	    << " is a face of " << simplex_format(cof, false) << "\n";
      ofs << "            </ul>\n";
      ofs << "            <ul>\n";
      if (cc_info.status_ == CC_detail_info::Result_type::self) {
	ofs << "            <p><span style=\"color:blue\">The simplex "
	    << simplex_format(cc_info.simplex_, false)
	    << " already exists in the cell complex!</span></p>\n";
      }
      if (cc_info.status_ == CC_detail_info::Result_type::face) {
	ofs << "            <p><span style=\"color:red\">The simplex "
	    << simplex_format(cc_info.simplex_, false) << " is a face of the simplex "
	    << simplex_format(cc_info.trigger_, false) << "!</span><br>\n";
	ofs << "              <ul>\n";
	for (const std::string post_face: cc_info.post_faces_)
	  ofs << "                <li id = \"" << id_from_simplex("I" + post_face) << "\">"
	      << "Post deleting " << simplex_format(post_face, false) << "</li>\n";
	ofs << "              </ul>\n";
	ofs << "            </p>\n";
	ofs << "            <p id = \"" << id_from_simplex("I" + cc_info.trigger_) << "\">"
	    << "Deleting " << simplex_format(cc_info.trigger_, false) << "</p>\n";
      }
      // for (const std::string& fac: cc_info.cofaces_)
      // 	ofs << "              <li>Checking if " << simplex_format(cc_info.simplex_, false)
      // 	    << " is a coface of " << simplex_format(fac, false) << "\n";
      if (cc_info.status_ == CC_detail_info::Result_type::coface) {
	ofs << "            <p><span style=\"color:darkorange\">The simplex "
	    << simplex_format(cc_info.simplex_, false) << " is a coface of the simplex "
	    << simplex_format(cc_info.trigger_, false) << "!</span><p>\n";
      }
      if (cc_info.status_ == CC_detail_info::Result_type::inserted) {
	ofs << "            <p><span style=\"color:green\">Successfully inserted "
	    << simplex_format(cc_info.simplex_, false) << "!</span><p>\n";
      }
      ofs << "            </ul>\n";
      ofs << "          </li>\n";
    }
    ofs << "        </ul>\n";

    if (i < cc_boundary_summary_lists.size()) {
      ofs << "      <h4 id=\"dim" << i << "b\"> Summary for boundary simplices</h4>\n";
      ofs << "      <p><a href=\"#dim" << i << "i\">Go to interior</a></p>\n";
      ofs << "        <ul>\n";
      for (const CC_summary_info& cc_info: cc_boundary_summary_lists[i])
	ofs << "          <li  id = \"" << id_from_simplex("B" + cc_info.face_) << "\">"
	    << simplex_format(cc_info.face_, true)
	    << " cell =" << cc_info.cell_ << "</li>\n";
      ofs << "        </ul>\n";
      ofs << "      <h4> Prejoin state of the boundary cells of dimension " << i << "</h4>\n";
      auto prejoin_it = cc_boundary_prejoin_lists[i].begin();
      while (prejoin_it != cc_boundary_prejoin_lists[i].end()) {
	std::size_t j = prejoin_it->dimension_;
	ofs << "        <h5>" << j << "-dimensional ambient simplices</h5>\n";
	ofs << "          <ul>\n";
	for (; prejoin_it->dimension_ == j; ++prejoin_it) {
	  ofs << "          <li>" << simplex_format(prejoin_it->simplex_, true)
	      << " join = " << simplex_format(prejoin_it->join_, true)
	      << " boundary:\n"
	      << "            <ul>\n";
	  for (const auto& face: prejoin_it->faces_)
	    ofs << "              <li>" << simplex_format(face) << "</li>";
	  ofs << "            </ul>\n";
	  switch (prejoin_it->status_) {
	  case (CC_prejoin_info::Result_type::join_single):
	    ofs << "            <p style=\"color: red\">Deleted "
		<< simplex_format(prejoin_it->simplex_, true)
		<< " as it has a single face.</p>";
	    break;
	  case (CC_prejoin_info::Result_type::join_is_face):
	    ofs << "            <p style=\"color: red\">Deleted "
		<< simplex_format(prejoin_it->simplex_, true)
		<< " as its join " << simplex_format(prejoin_it->join_, true)
		<< " is one of the faces.</p>";
	    break;
	  case (CC_prejoin_info::Result_type::join_different):
	    ofs << "            <p style=\"color: magenta\">Deleted " << simplex_format(prejoin_it->simplex_, true)
		<< " and replaced by its join " << simplex_format(prejoin_it->join_, true)
		<< ".</p>";
	    break;
	  case (CC_prejoin_info::Result_type::join_same):
	    ofs << "            <p style=\"color: green\">Kept " << simplex_format(prejoin_it->simplex_, true)
		<< ".</p>";
	  }
	  ofs << "          </li>";
	}
	ofs << "          </ul>\n";
      }
    }
    if (i < cc_boundary_detail_lists.size()) {
      ofs << "      <h4> Details for boundary simplices</h4>\n"
	  << "        <ul>\n";
      for (const CC_detail_info& cc_info: cc_boundary_detail_lists[i]) {
	if (cc_info.status_ == CC_detail_info::Result_type::join_single) {
	  ofs << "          <li style=\"color:magenta\" id = \""
	      << id_from_simplex("B" + cc_info.simplex_) << "\"> Simplex "
	      << simplex_format(cc_info.simplex_, true) << " has only one face ("
	      << simplex_format(cc_info.trigger_, true) << ") and is deleted.";
	  continue;
	}
	if (cc_info.status_ == CC_detail_info::Result_type::join_single) {
	  ofs << "          <li style=\"color:darkmagenta\" id = \""
	      << id_from_simplex("B" + cc_info.simplex_) << "\"> The join of the simplex "
	      << simplex_format(cc_info.simplex_, true) << " is one of its faces ("
	      << simplex_format(cc_info.trigger_, true) << "), hence it is is deleted.";
	  continue;
	}
	ofs << "          <li> Insert_simplex called on " << simplex_format(cc_info.simplex_, true);
	ofs << "            <ul>\n";
	// for (const std::string& cof: cc_info.faces_)
	//   ofs << "              <li>Checking if " << simplex_format(cc_info.simplex_, true)
	//       << " is a face of " << simplex_format(cof, true) << "\n";
	ofs << "            </ul>\n";
	ofs << "            <ul>\n";
	if (cc_info.status_ == CC_detail_info::Result_type::self) {
	  ofs << "            <p><span style=\"color:blue\">The simplex "
	      << simplex_format(cc_info.simplex_, true)
	      << " already exists in the cell complex!</span></p>\n";
	}
	if (cc_info.status_ == CC_detail_info::Result_type::face) {
	  ofs << "            <p><span style=\"color:red\">The simplex "
	      << simplex_format(cc_info.simplex_, true) << " is a face of the simplex "
	      << simplex_format(cc_info.trigger_, true) << "!</span><br>\n";
	  ofs << "              <ul>\n";
	  for (const std::string post_face: cc_info.post_faces_)
	    ofs << "                <li id=\"" << id_from_simplex("B" + post_face)
		<< "\">Post deleting " << simplex_format(post_face, true)
		<< "</li>\n";
	  ofs << "              </ul>\n";
	  ofs << "            </p>\n";
	  ofs << "            <p id=\"" << id_from_simplex(cc_info.trigger_)
	      << "\">Deleting " << simplex_format(cc_info.trigger_, true) << "</p>\n";
	}
	// for (const std::string& fac: cc_info.cofaces_)
	//   ofs << "              <li>Checking if " << simplex_format(cc_info.simplex_, true)
	//       << " is a coface of " << simplex_format(fac, true) << "\n";
	ofs << "            </ul>\n";
	ofs << "          </li>\n";
	if (cc_info.status_ == CC_detail_info::Result_type::coface) {
	  ofs << "            <p><span style=\"color:darkorange\">The simplex "
	      << simplex_format(cc_info.simplex_, true) << " is a coface of the simplex "
	      << simplex_format(cc_info.trigger_, true) << "!</span><p>\n";
	}
	if (cc_info.status_ == CC_detail_info::Result_type::inserted) {
	  ofs << "            <p><span style=\"color:green\">Successfully inserted "
	      << simplex_format(cc_info.simplex_, true) << "!</span><p>\n";
	}
      }
      ofs << "        </ul>\n";
    }
  }  
  ofs << "    </div>\n";
}

void write_visu(std::ofstream& ofs) {
  ofs << "    <div id=\"visu\">\n"
      << "      <h2> Visualization details debug trace </h2>\n";
  // std::vector<std::map<std::string, std::string> > vs_maps(cc_interior_summary_lists.size());
  std::map<std::string, std::string> vs_map;
  for (const auto& sv_pair: simplex_vlist_map)
    vs_map.emplace(std::make_pair(sv_pair.second, sv_pair.first));
  ofs << "      <ul>\n";
  for (const auto& vs_pair: vs_map) {
    std::string w_simplex = vs_pair.second.substr(1);
    bool is_boundary = vs_pair.second[0] == 'B';
    ofs << "      <li><b>" << vs_pair.first << "</b>: "
	<< simplex_format(w_simplex, is_boundary) << "</li>\n";
  }
  ofs << "      </ul>\n";
  ofs << "    </div>\n";
}

void write_to_html(std::string file_name) {
  std::ofstream ofs(file_name + ".html", std::ofstream::out);
  ofs << "<!DOCTYPE html>\n"
      << "<html>\n";
  write_head(ofs);
  ofs << "  <body>\n";
  write_nav(ofs);
  ofs << "    <h1> Debug traces for " << file_name << " </h1>\n";
  write_mt(ofs);
  write_cc(ofs);
  write_visu(ofs);  
  ofs << "  </body>\n";
  ofs << "</html>\n";

  ofs.close();
}

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif

