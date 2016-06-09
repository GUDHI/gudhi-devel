#ifndef OUTPUT_TIKZ_H
#define OUTPUT_TIKZ_H

#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <cmath>

typedef double FT;


//////////////////AUX/////////////////////

void write_preamble(std::ofstream& ofs)
{
  ofs << "\\documentclass{standalone}\n"
      << "\\usepackage{tikz}\n\n"
      << "\\begin{document}\n"
      << "\\begin{tikzpicture}\n";
}

void write_end(std::ofstream& ofs)
{
  ofs << "\\end{tikzpicture}\n"
      << "\\end{document}";
}


/////////////////MAIN//////////////////////

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


// A little script to make a tikz histogram of epsilon distribution
// Returns the average epsilon
void write_histogram(std::vector<double> histo, std::string file_name = "histogram.tikz", std::string xaxis = "$\\epsilon/\\epsilon_{max}$",  std::string yaxis = "$\\epsilon$", FT max_x = 1)
{
  int n = histo.size();
  
  std::ofstream ofs (file_name, std::ofstream::out);
  FT barwidth = 20.0/n;
  FT max_value = *(std::max_element(histo.begin(), histo.end()));
  std::cout << max_value << std::endl;
  FT ten_power = pow(10, ceil(log10(max_value)));
  FT max_histo = ten_power;
  if (max_value/ten_power > 1) {
    if (max_value/ten_power < 2)
      max_histo = 0.2*ten_power;
    else if (max_value/ten_power < 5)
      max_histo = 0.5*ten_power;
  }
  std::cout << ceil(log10(max_value)) << std::endl << max_histo << std::endl;
  FT unitht = max_histo/10.0;
  write_preamble(ofs);
  
  ofs << "\\draw[->] (0,0) -- (0,11);\n" <<
    "\\draw[->] (0,0) -- (21,0);\n" <<
    "\\foreach \\i in {1,...,10}\n" <<
    "\\draw (0,\\i) -- (-0.1,\\i);\n" <<
    "\\foreach \\i in {1,...,20}\n" <<
    "\\draw (\\i,0) -- (\\i,-0.1);\n" <<
 
    "\\node at (-1,11) {" << yaxis << "};\n" << 
    "\\node at (22,-1) {" << xaxis << "};\n" << 
    "\\node at (-0.5,-0.5) {0};\n" << 
    "\\node at (-0.5,10) {" << max_histo << "};\n" << 
    "\\node at (20,-0.5) {" << max_x << "};\n";
    
  for (int i = 0; i < n; ++i)
    ofs << "\\draw (" << barwidth*i << "," << histo[i]/unitht << ") -- ("
        << barwidth*(i+1) << "," << histo[i]/unitht << ") -- ("
  << barwidth*(i+1) << ",0) -- (" << barwidth*i << ",0) -- cycle;\n";

  write_end(ofs);  
  ofs.close();
}

struct Pers_interval {
  double alpha_start, alpha_end;
  int dim;
  Pers_interval(double alpha_start_, double alpha_end_, int dim_)
    : alpha_start(alpha_start_), alpha_end(alpha_end_), dim(dim_)
  {}
};

void write_barcodes(std::string in_file, double alpha2, std::string out_file = "barcodes.tikz.tex")
{
  std::ifstream ifs(in_file, std::ios::in);
  std::string line;
  std::vector<Pers_interval> pers_intervals;
  while (getline(ifs, line)) {
    int p, dim;
    double alpha_start, alpha_end;
    std::istringstream iss(line);
    iss >> p >> dim >> alpha_start >> alpha_end;
    if (alpha_start != alpha_end) {
      if (alpha_end < alpha_start)
        alpha_end = alpha2;
      pers_intervals.push_back(Pers_interval(alpha_start, alpha_end, dim));
    }
  }
  ifs.close();
  std::ofstream ofs (out_file, std::ofstream::out);
  write_preamble(ofs);
  double barwidth = 0.01;
  int i = 0;
  for (auto interval: pers_intervals) {
    std::string color = "black";
    switch (interval.dim) {
    case 0: color = "orange"; break;
    case 1: color = "red"; break;
    case 2: color = "blue"; break;
    case 3: color = "green"; break;
    case 4: color = "yellow"; break;
    default: color = "orange"; break;
    }
    ofs << "\\fill[" << color << "] (" << interval.alpha_start << "," << barwidth*i << ") rectangle ("
        << interval.alpha_end << "," << barwidth*(i+1) <<");\n";
    i++;
  }
  write_end(ofs);
  ofs.close();
}

#endif
