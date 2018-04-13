#ifndef OUTPUT_POINTS_TO_MEDIT_
#define OUTPUT_POINTS_TO_MEDIT_

template <class Point_range>          
void output_points_to_medit(Point_range& range, std::string file_name = "points.mesh") {
  short d = range.begin()->size();
  if (d > 3);

  std::ofstream ofs (file_name, std::ofstream::out);
  if (d <= 2)
    ofs << "MeshVersionFormatted 1\nDimension 2\n";
  else
    ofs << "MeshVersionFormatted 1\nDimension 3\n";

  ofs << "Vertices\n" << range.size() << "\n";
  for (auto v: range) {
    if (d == 2)
      ofs << v[0] << " " << v[1] << " ";
    else
      ofs << v[0] << " " << v[1] << " " << v[2] << " ";
    ofs << "512\n";
  }
  
  // ofs << "Vertices\n" << (d+1)*range.size() << "\n";
  // for (auto v: range) 
  //   for (int i = 0; i < d+1; ++i) {
  //     for (auto x: v)
  //       ofs << x << " ";
  //     ofs << "\n";
  //   }
  if (d == 2) {
    // ofs << "Edges\n" << 3*range.size() << "\n";
    // for (unsigned k = 0; k < range.size(); ++k) {
    //   ofs << 3*k+1 << " " << 3*k+2 << "\n";
    //   ofs << 3*k+1 << " " << 3*k+3 << "\n";
    //   ofs << 3*k+2 << " " << 3*k+3 << "\n";
    // }
    // ofs << "Triangles\n" << range.size() << "\n";
    // for (unsigned k = 0; k < range.size(); ++k)
    //   ofs << 3*k+1 << " " << 3*k+2 << " " << 3*k+3 << "\n";
  }

  ofs.close();
}

#endif
