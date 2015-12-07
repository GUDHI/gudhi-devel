#ifndef OUTPUT_TIKZ_H
#define OUTPUT_TIKZ_H

#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <cmath>

void write_tikz_plot(std::vector<FT> data, std::string filename)
{
  int n = data.size();
  FT vmax = *(std::max_element(data.begin(), data.end()));
  //std::cout << std::log10(vmax) << " " << std::floor(std::log10(vmax));
  
  FT order10 = pow(10,std::floor(std::log10(vmax)));
  int digit = std::floor( vmax / order10) + 1;
  if (digit == 4 || digit == 6) digit = 5;
  if (digit > 6) digit = 10;
  FT plot_max = digit*order10;
  std::cout << plot_max << " " << vmax;
  FT hstep = 10.0/(n-1);
  FT wstep = 10.0 / plot_max;

  std::cout << "(eps_max-eps_min)/(N-48) = " << (vmax-*data.begin())/(data.size()-48) << "\n";
  std::ofstream ofs(filename, std::ofstream::out);

  ofs <<
    "\\documentclass{standalone}\n" <<
    "\\usepackage[utf8]{inputenc}\n" <<
    "\\usepackage{amsmath}\n" <<
    "\\usepackage{tikz}\n\n" <<
    "\\begin{document}\n" <<
    "\\begin{tikzpicture}\n";
  
  ofs <<  "\\draw[->] (0,0) -- (0,11);" << std::endl <<
    "\\draw[->] (0,0) -- (11,0);" << std::endl <<
    "\\foreach \\i in {1,...,10}" << std::endl <<
    "\\draw (0,\\i) -- (-0.05,\\i);" << std::endl <<
    "\\foreach \\i in {1,...,10}" << std::endl <<
    "\\draw (\\i,0) -- (\\i,-0.05);" << std::endl << std::endl <<

    "\\foreach \\i in {1,...,10}" << std::endl <<
    "\\draw[dashed] (-0.05,\\i) -- (11,\\i);" << std::endl << std::endl <<
    
    "\\node at (-0.5,11) {$*$}; " << std::endl <<
    "\\node at (11,-0.5) {$*$}; " << std::endl <<
    "\\node at (-0.5,-0.5) {0}; " << std::endl <<
    "\\node at (-0.5,10) {" << plot_max << "}; " << std::endl <<
    "%\\node at (10,-0.5) {2}; " << std::endl;

  ofs << "\\draw[red] (0," << wstep*data[0] << ")";
  for (int i = 1; i < n; ++i)
    ofs << " -- (" << hstep*i << "," << wstep*data[i] << ")";
  ofs << ";\n";

  ofs <<
    "\\end{tikzpicture}\n" <<
    "\\end{document}";
  
  ofs.close();


  
}

#endif
