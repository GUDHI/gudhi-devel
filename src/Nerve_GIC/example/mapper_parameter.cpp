#include "mapper.h"
#define ETA 0.001
#define VERBOSE 1
#define DELTA 1
#define RESOLUTION 1
#define CONSTANT 10

double GetUniform(){
  static default_random_engine re;
  static uniform_real_distribution<double> Dist(0,1);
  return Dist(re);
}

void SampleWithoutReplacement(int populationSize, int sampleSize, vector<int> & samples){
  int& n = sampleSize; int& N = populationSize;
  int t = 0; int m = 0; double u;
  while (m < n){
    u = GetUniform();
    if ( (N - t)*u >= n - m )
      t++;
    else{samples[m] = t; t++; m++;}
  }
}

int main(int argc, char** argv){

    if(argc <= 7 || argc >= 9){cout << "./param <cloud_file> <function:/coordinate:/eccentricity:> <func_file/number/matrix_file>" << \
                                        " <VNE> <gain> <subsampling:/graph:> <N/graph_name>" << endl; return 0;}

    char* const cloud = argv[1];
    char* const funct = argv[2];
    bool VNE = atoi(argv[4]);
    double g = atof(argv[5]);
    char* const method = argv[6];

    Cloud C;

    if(VERBOSE)  cout << "Reading input cloud from file " << cloud << "..." << endl;
    C = read_cloud(cloud); int n = C.size();
    if(VERBOSE)  cout << "  Done." << endl;

    if (strcmp(funct, "function:") == 0){
      char* const func = argv[3];
      if(VERBOSE)  cout << "Reading input filter from file " << func << "..." << endl;
      read_function_from_file(func,C);
    }
    if (strcmp(funct, "coordinate:") == 0){
      int number = atoi(argv[3]);
      if(VERBOSE)  cout << "Using coordinate " << number << " as a filter..." << endl;
      read_coordinate(number,C);
    }
    if (strcmp(funct, "eccentricity:") == 0){
      char* const matrix = argv[3];
      cout << "Computing eccentricity with distance matrix " << matrix << "..." << endl;
      compute_eccentricity(C,matrix);
    }
    if(VERBOSE)  cout << "  Done." << endl;

    double maxf = C[0].func; double minf = C[0].func;
    for(int i = 1; i < n; i++){ maxf = max(maxf,C[i].func); minf = min(minf,C[i].func);}

    double delta = 0; double lip = 0;

    if (strcmp(method, "subsampling:") == 0){

      char name[100];
      if(VNE)  sprintf(name,"%s_VNdist",cloud);
      else  sprintf(name,"%s_dist",cloud);
      ifstream input(name, ios::out | ios::binary);
      vector<vector<double> > dist; dist.clear(); double d;

      if(input.good()){
        cout << "Reading distances..." << endl;
        for(int i = 0; i < n; i++){
          vector<double> dis; dis.clear();
          for (int j = 0; j < n; j++){input.read((char*) &d,8); dis.push_back(d);}
          dist.push_back(dis);
        }
        input.close();
      }
      else{
        cout << "Computing distances..." << endl; vector<double> V;
        input.close(); ofstream output(name, ios::out | ios::binary);
        if(VNE){
          int dim = C[0].coord.size();
          for(int i = 0; i < dim; i++){
            double meani = 0;
            for (int j = 0; j < n; j++)  meani += C[j].coord[i]/n;
            double vari = 0;
            for (int j = 0; j < n; j++)  vari += pow(C[j].coord[i]-meani,2)/n;
            V.push_back(vari);
          }
        }

        for(int i = 0; i < n; i++){
          if( (int) floor( 100*((double) i)/((double) n)+1 ) %10 == 0  ){cout << "\r" << floor( 100*((double) i)/((double) n) +1) << "%" << flush;}
          vector<double> dis; dis.clear();
          for (int j = 0; j < n; j++){
            if(!VNE){d = C[i].EuclideanDistance(C[j]); dis.push_back(d);}
            else{d = C[i].VarianceNormalizedEuclideanDistance(C[j],V); dis.push_back(d);}
            output.write((char*) &d,8);
          }
          dist.push_back(dis);
        }
        output.close();
      }

      cout << "Done." << endl;

      int m = floor(n/pow(log(n)/log(CONSTANT),1+ETA)); m = min(m,n-1);

      if(VERBOSE)  cout << "dimension = " << C[0].coord.size() << endl << "n = " << n << " and s(n) = " << m << endl;
      if(VERBOSE)  cout << "range = [" << minf << ", " << maxf << "]" << endl;

      vector<int> samples(m); int N = atoi(argv[7]);
      #pragma omp parallel for
      for (int i = 0; i < N; i++){
        SampleWithoutReplacement(n,m,samples);
        double hausdorff_dist = 0;
        for (int j = 0; j < n; j++){
          Point p = C[j];
          double mj = dist[j][samples[0]];
          double Mi = abs(p.func-C[0].func)/(dist[j][0]);
          for (int k = 1; k < m; k++){mj = min(mj, dist[j][samples[k]]);}
          for (int k = j+1; k < n; k++){Mi = max(Mi, abs(p.func-C[k].func)/(dist[j][k]));}
          hausdorff_dist = max(hausdorff_dist, mj); lip = max(lip,Mi);
        }
        delta += hausdorff_dist;
      }
      delta /= N;

    }

    if (strcmp(method, "graph:") == 0){
      char* const graph_name = argv[7];
      cout << "Reading neighborhood graph from file " << graph_name << "..." << endl;
      pair<double,double> P = compute_delta_from_file(C,graph_name);
      delta = P.first; lip = P.second;
    }

    if(VERBOSE)  cout << "lip = " << lip << endl;

    if(VERBOSE)  cout << "delta = " << delta << endl;
    else cout << delta << endl;

    if(VERBOSE)  cout << "g = " << g << endl;
    else  cout << g << endl;

    if(!DELTA){
      if(VERBOSE){
        if(RESOLUTION)  cout << "r = " << lip*delta/g << endl;
        else  cout << "r = " << floor(g*(maxf-minf)/(lip*delta*(1-g))) << endl;}
      else{
        if(RESOLUTION)  cout << lip*delta/g << endl;
        else  cout << floor(g*(maxf-minf)/(lip*delta*(1-g))) << endl;}}
    else{
      if(VERBOSE){
        if(RESOLUTION)  cout << "r = " << lip*delta << endl;
        else  cout << "r = " << floor((maxf-minf)/(lip*delta*(1-g))) << endl;}
      else{
        if(RESOLUTION)  cout << lip*delta << endl;
        else  cout << floor((maxf-minf)/(lip*delta*(1-g))) << endl;}}

    return 0;
}
