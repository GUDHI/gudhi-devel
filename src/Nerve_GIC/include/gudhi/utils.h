#include <cstdlib>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <limits>
#include <assert.h>
#include <math.h>
#include <iterator>
#include <omp.h>
#include <list>
#include <random>

using namespace std;

class Point{

  public:

    int ID;
    double func;
    vector<double> coord;

    Point(){ ID = -1;}
    Point(const int & _ID){ID = _ID;} // constructor 1
    Point(const int & _ID, const vector<double> & _coord){ID = _ID; coord = _coord;} // constructor 2
    bool operator<(const Point & p) const {if(func != p.func){return func < p.func;}else{return ID < p.ID;}} // comparator
    bool operator==(const Point & p) const {return ID == p.ID;}
    double EuclideanDistance(const Point & p) const {
      //cout << this->coord.size() << " " << p.coord.size() << endl;
      assert (coord.size() == p.coord.size());
      double x = 0; int dim = p.coord.size(); for (int i = 0; i < dim; i++){x+=(coord[i]-p.coord[i])*(coord[i]-p.coord[i]);}
      return sqrt(x);
    }
    double VarianceNormalizedEuclideanDistance(const Point & p, const vector<double> & V) const {
      //cout << this->coord.size() << " " << p.coord.size() << endl;
      assert (coord.size() == p.coord.size());
      double x = 0; int dim = p.coord.size(); for (int i = 0; i < dim; i++){x+=(coord[i]-p.coord[i])*(coord[i]-p.coord[i])/V[i];}
      return sqrt(x);
    }

};

typedef vector<Point> Cloud;

Cloud read_cloud(char* const & name){
  int dim; Cloud C; C.clear();
  ifstream input(name);
  if(input){
    string line; vector<double> coord; int ID = 0;
    while(getline(input,line)){
      coord.clear(); double x = numeric_limits<double>::max();
      stringstream stream(line);
      while(stream >> x){coord.push_back(x);}
      Point p = Point(ID, coord);
      if(x != numeric_limits<double>::max()){C.push_back(p); ID++;}
    }
  }
  else{cout << "Failed to read file " << name << endl; return C;}
  return C;
}

void read_function_from_file(char* const & name, Cloud & C){
  int num_pts = C.size(); double x;
  ifstream input(name);
  if(input){
    for(int i = 0; i < num_pts; i++)
      input >> C[i].func;
  }
  else{cout << "Failed to read file " << name << endl; return;}
  return;
}

void read_coordinate(const int & number, Cloud & C){
  int num_pts = C.size();
  for(int i = 0; i < num_pts; i++)
    C[i].func = C[i].coord[number];
  return;
}


