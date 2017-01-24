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

double EuclideanDistance(const vector<double> & V1, const vector<double> & V2){
  assert(V1.size() == V2.size());
  double res = 0;
  for(int i = 0; i < V1.size(); i++)  res += pow(V1[i]-V2[i],2);
  return sqrt(res);
}

int main(int argc, char** argv){

  char* const cloud = argv[1];
  char* const method = argv[2];

  if(argc <= 4){cout << "Usage: <cloud_file> <DTM/GDE> <k/h> <threshold>" << endl; return 0;}
  if(argc >= 6){cout << "Too many arguments !" << endl; return 0;}

  ifstream Cloud(cloud);
  vector<vector<double> > C, distances; C.clear(); double d;
  string line;

  int nb = 0; cout << "Reading cloud..." << endl;
  while(getline(Cloud,line)){
    if(nb%100000 == 0)  cout << "  " << nb << endl;
    vector<double> P;
    stringstream stream(line);
    while(stream >> d)  P.push_back(d);
    C.push_back(P); nb++;
  }
  cout << endl;

  for(int i = 0; i < nb; i++){vector<double> vect(nb); distances.push_back(vect);}

  char name[100]; sprintf(name,"%s_dist",cloud);
  ifstream input(name, ios::out | ios::binary);

  if(input.good()){

    cout << "  Reading distances..." << endl;

    for(int i = 0; i < nb; i++){
      for (int j = 0; j < nb; j++){
        input.read((char*) &d,8); distances[i][j] = d;
      }
    }

    input.close();
    cout << "  Done." << endl;

  }

  else{

    cout << "  Computing distances..." << endl;
    input.close(); ofstream output(name, ios::out | ios::binary);

    for(int i = 0; i < nb-1; i++){
      if( (int) floor( 100*((double) i)/((double) nb)+1 ) %10 == 0  ){cout << "\r" << floor( 100*((double) i)/((double) nb) +1) << "%" << flush;}
      for (int j = i+1; j < nb; j++){
        d = EuclideanDistance(C[i],C[j]); distances[i][j] = d; distances[j][i] = d;
        output.write((char*) &d,8);
      }
    }

    output.close();
    cout << endl << "  Done." << endl;

  }

  if(strcmp(method,"DTM") == 0){

    int k = atoi(argv[3]);
    double thresh = atof(argv[4]);

    char nameo[100]; sprintf(nameo,"%s_DTM_%d_%g",cloud,k,thresh);
    ofstream output(nameo); vector<double> densities; double M = 0;

    cout << "Computing DTM..." << endl;
    for(int i = 0; i < nb; i++){
      if( (int) floor( 100*((double) i)/((double) nb)+1 ) %10 == 0  ){cout << "\r" << floor( 100*((double) i)/((double) nb) +1) << "%" << flush;}
      double density = 0;
      sort(distances[i].begin(),distances[i].end());
      for(int j = 0; j < k; j++)
        density += pow(distances[i][j],2);
      density = sqrt(k)/sqrt(density);
      if(density >= M)  M = density;
      densities.push_back(density);

    }
    cout << endl << "  Done." << endl;

    for(int i = 0; i < nb; i++){
      if(densities[i] >= thresh*M){
        for(int j = 0; j < C[i].size(); j++)  output << C[i][j] << " ";
        output << endl;
      }
    }

    output.close();

  }

  if(strcmp(method,"GaussianDE") == 0){

    double bandwidth = atof(argv[3]);
    double thresh = atof(argv[4]);

    char nameo[100]; sprintf(nameo,"%s_GDE_%g_%g",cloud,bandwidth,thresh);
    ofstream output(nameo); vector<double> densities; double M = 0;

    cout << "Computing GDE..." << endl;
    for(int i = 0; i < nb; i++){
      if( (int) floor( 100*((double) i)/((double) nb)+1 ) %10 == 0  ){cout << "\r" << floor( 100*((double) i)/((double) nb) +1) << "%" << flush;}
      double density = 0;
      for(int j = 0; j < nb; j++)
        density += exp(-pow(distances[i][j],2)/(2*pow(bandwidth,2)));
      density /= sqrt(2*3.1415)*nb*bandwidth;
      if(density >= M)  M = density;
      densities.push_back(density);

    }
    cout << endl << "  Done." << endl;

    for(int i = 0; i < nb; i++){
      if(densities[i] >= thresh*M){
        for(int j = 0; j < C[i].size(); j++)  output << C[i][j] << " ";
        output << endl;
      }
    }

    output.close();

  }

  return 0;

}
