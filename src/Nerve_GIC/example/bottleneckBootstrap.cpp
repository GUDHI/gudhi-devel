 #include "mapper.h"
#define ETA 0.001
#define SHIFT 1
#define DELTA 0
#define CONSTANT 10

double GetUniform(){
  static default_random_engine re;
  static uniform_real_distribution<double> Dist(0,1);
  return Dist(re);
}

void SampleWithReplacement(int populationSize, int sampleSize, vector<int> & samples){
  int m = 0; double u;
  while (m < sampleSize){
    u = GetUniform();
    samples[m] = floor(u*populationSize);
    m++;
  }
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

void compute_estimator(Cloud & C, char* const & cloud, const bool & VNE, const int & Nsubs, \
                       const double & gain){

    int n = C.size();
    double** dist = new double*[n];
    for(int i = 0; i < n; i++)  dist[i] = new double[n];
    double d;
    Cover I; AdjacencyMatrix G; MapperElements M;
    double maxf = C[0].func; double minf = C[0].func;
    for(int i = 1; i < n; i++){ maxf = max(maxf,C[i].func); minf = min(minf,C[i].func);}

    char name1[100];
    if(VNE)  sprintf(name1, "%s_VNdist", cloud);
    else  sprintf(name1, "%s_dist", cloud);
    char* const name2 = name1;
    ifstream input(name2, ios::out | ios::binary);

    if(input.good()){
      cout << "  Reading distances for parameter selection..." << endl;
      for(int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
          input.read((char*) &d,8); dist[i][j] = d;
        }
      }
      input.close();
    }
    else{
      cout << "  Computing distances for parameter selection..." << endl; vector<double> V;
      input.close(); ofstream output(name2, ios::out | ios::binary);
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
        if( (int) floor( 100*((double) i)/((double) n)+1 ) %10 == 0  ){
          cout << "\r" << floor( 100*((double) i)/((double) n) +1) << "%" << flush;}
        for (int j = 0; j < n; j++){
          if(!VNE){d = C[i].EuclideanDistance(C[j]); dist[i][j] = d;}
          else{d = C[i].VarianceNormalizedEuclideanDistance(C[j],V); dist[i][j] = d;}
          output.write((char*) &d,8);
        }
      }
      output.close();
    }

    int m = floor(n/pow(log(n)/log(CONSTANT),1+ETA));
    double lip = 0; double delta = 0; vector<int> samples(m);
    #pragma omp parallel for
    for (int i = 0; i < Nsubs; i++){
      SampleWithoutReplacement(n,m,samples);
      double hausdorff_dist = 0;
      for (int j = 0; j < n; j++){
        Point p = C[j]; double mj = dist[j][samples[0]]; double Mi = abs(p.func-C[0].func)/dist[j][0];
        for (int k = 1; k < m; k++)  mj = min(mj, dist[j][samples[k]]);
        if(j < n-1)
          for (int k = j+1; k < n; k++)  Mi = max(Mi, abs(p.func-C[k].func)/dist[j][k]);
        hausdorff_dist = max(hausdorff_dist, mj); lip = max(lip,Mi);
      }
      delta += hausdorff_dist;
    }
    delta /= Nsubs; double resolution; if(DELTA)  resolution = lip*delta; else  resolution = lip*delta/gain;
    //cout << delta << " " << lip << " " << resolution << endl;
    cout << "  Done." << endl;

    G = build_neighborhood_graph_from_scratch_with_parameter(C,delta,name2,VNE);

    sort(C.begin(),C.end());
    I = Cover(C.begin()->func, (C.end()-1)->func, resolution, gain, SHIFT, 1);
    I.proper_value();

    map<int,double> colors; colors.clear(); map<int,int> num; num.clear();
    map<Point,vector<int> > clusters; vector<int> dum; dum.clear();
    for(int i = 0; i < C.size(); i++)  clusters.insert(pair<Point,vector<int> >(C[i],dum)); vector<double> col(G.size());

    M = MapperElts(G,I,clusters,col);
    SparseGraph MG;
    if (DELTA)  MG = build_mapperDelta_graph(M,G,clusters,colors,num); else  MG = build_mapper_graph(M,colors,num);

    char mappg[100]; sprintf(mappg,"%s_mapper.txt",cloud);
    ofstream graphicg(mappg);

    double maxv, minv; maxv = numeric_limits<double>::min(); minv = numeric_limits<double>::max();
    for (map<int,double>::iterator iit = colors.begin(); iit != colors.end(); iit++){
      if(iit->second > maxv)  maxv = iit->second;
      if(minv > iit->second)  minv = iit->second;
    }

    int k = 0; int ke = 0;

    for (SparseGraph::iterator it = MG.begin(); it != MG.end(); it++){
      k++; for (int j = 0; j < it->second.size(); j++)  ke++;
    }

    graphicg << k << " " << ke << endl; int kk = 0;
    for (SparseGraph::iterator it = MG.begin(); it != MG.end(); it++){graphicg << kk << " " << colors[it->first] << endl; kk++;}
    for (SparseGraph::iterator it = MG.begin(); it != MG.end(); it++)
      for (int j = 0; j < it->second.size(); j++)
        graphicg << it->first << " " << it->second[j] << endl;
}

void compute_bootstrapped_estimator(const vector<int> & samp, char* const & cloud, const int & Nboot, double** dist, Cloud & Cboot, \
                                           const int & Nsubs, const double & gain){

    int n = Cboot.size();
    Cover I; MapperElements M;
    double maxf = Cboot[0].func; double minf = Cboot[0].func;
    for(int i = 1; i < n; i++){ maxf = max(maxf,Cboot[i].func); minf = min(minf,Cboot[i].func);}
    int m = floor(n/pow(log(n)/log(CONSTANT),1+ETA));
    double lip = 0; double delta = 0; vector<int> samples(m);
    #pragma omp parallel for
    for (int i = 0; i < Nsubs; i++){
      SampleWithoutReplacement(n,m,samples);
      double hausdorff_dist = 0;
      for (int j = 0; j < n; j++){
        Point p = Cboot[j]; double mj = dist[samp[j]][samp[samples[0]]]; double Mi = abs(p.func-Cboot[0].func)/dist[samp[j]][samp[0]];
        for (int k = 1; k < m; k++)  mj = min(mj, dist[samp[j]][samp[samples[k]]]);
        if(j < n-1)
          for (int k = j+1; k < n; k++)  Mi = max(Mi, abs(p.func-Cboot[k].func)/dist[samp[j]][samp[k]]);
        hausdorff_dist = max(hausdorff_dist, mj); lip = max(lip,Mi);
      }
      delta += hausdorff_dist;
    }
    delta /= Nsubs; double resolution;
    if(DELTA)  resolution = lip*delta; else  resolution = lip*delta/gain;
    //cout << delta << " " << lip << " " << resolution << endl;

    AdjacencyMatrix G; G.clear(); vector<Point> adj;

    for(int i = 0; i < n; i++){
      adj.clear();
      for(int j = 0; j < n; j++)
        if(dist[samp[i]][samp[j]] <= delta)
          adj.push_back(Cboot[j]);
      G.insert(pair<Point, vector<Point> >(Cboot[i],adj));
    }

    sort(Cboot.begin(),Cboot.end());
    I = Cover(Cboot.begin()->func, (Cboot.end()-1)->func, resolution, gain, SHIFT, 1);
    I.proper_value();

    map<int,double> colors; colors.clear(); map<int,int> num; num.clear();
    map<Point,vector<int> > clusters; vector<int> dum; dum.clear();
    for(int i = 0; i < Cboot.size(); i++)  clusters.insert(pair<Point,vector<int> >(Cboot[i],dum)); vector<double> col(n);

    M = MapperElts(G,I,clusters,col);
    SparseGraph MG;
    if(DELTA)  MG = build_mapperDelta_graph(M,G,clusters,colors,num); else  MG = build_mapper_graph(M,colors,num);

    char mappg[100]; sprintf(mappg,"%s_mapper_%d.txt",cloud,Nboot);
    ofstream graphicg(mappg);

    double maxv, minv; maxv = numeric_limits<double>::min(); minv = numeric_limits<double>::max();
    for (map<int,double>::iterator iit = colors.begin(); iit != colors.end(); iit++){
      if(iit->second > maxv)  maxv = iit->second;
      if(minv > iit->second)  minv = iit->second;
    }

    int k = 0; int ke = 0;

    for (SparseGraph::iterator it = MG.begin(); it != MG.end(); it++){
      k++; for (int j = 0; j < it->second.size(); j++)  ke++;
    }

    graphicg << k << " " << ke << endl; int kk = 0;
    for (SparseGraph::iterator it = MG.begin(); it != MG.end(); it++){graphicg << kk << " " << colors[it->first] << endl; kk++;}
    for (SparseGraph::iterator it = MG.begin(); it != MG.end(); it++)
      for (int j = 0; j < it->second.size(); j++)
        graphicg << it->first << " " << it->second[j] << endl;
    graphicg.close();

}


int main(int argc, char** argv){

    if(argc <= 7 || argc >= 9){
      cout << "./bottleBoot <cloud_file> <VNE> <function:/coordinate:/eccentricity:> <func_file/number/distance matrix> <gain> <Nsubs> <Nboot>" << endl;
      return 0;}

    char* const cloud = argv[1]; bool normalized = atoi(argv[2]); char* const funct = argv[3]; double gain = atof(argv[5]);
    int Nsubs = atoi(argv[6]); int Nboot = atoi(argv[7]);

    Cloud C, D, Cboot; C = read_cloud(cloud); int n = C.size();

    if (strcmp(funct, "function:") == 0){
      char* const func = argv[4]; read_function_from_file(func,C);
    }
    if (strcmp(funct, "coordinate:") == 0){
      int number = atoi(argv[4]); read_coordinate(number,C);
    }
    if (strcmp(funct, "eccentricity") == 0){
      char* const matrix = argv[4]; compute_eccentricity(C,matrix);
    }

    D = C;

    cout << "Computing estimator..." << endl;
    compute_estimator(C,cloud,normalized,Nsubs,gain);
    cout << "Estimator computed." << endl;

    char nam[100]; sprintf(nam,"%s_mapper.txt",cloud); char dg[100]; sprintf(dg,"%s_ExDg",cloud);
    char command[100];
    sprintf(command,"cat %s | ../persistence/filtgraph 1 | ../persistence/pers 1 >> %s",nam,dg); system(command);
    char list[100]; sprintf(list,"%s_quantile",cloud);

    char name[100];
    if(normalized)  sprintf(name, "%s_VNdist", cloud);
    else  sprintf(name, "%s_dist", cloud);
    ifstream input(name, ios::out | ios::binary); double d;
    double** dist = new double*[n]; for(int i = 0; i < n; i++)  dist[i] = new double[n];
    for(int i = 0; i < n; i++)  for (int j = 0; j < n; j++){input.read((char*) &d,8); dist[i][j] = d;}
    input.close();

    for(int i = 0; i < Nboot; i++){

      Cboot.clear();
      vector<int> samples(n); SampleWithReplacement(n,n,samples); sort(samples.begin(), samples.end());
      for (int i = 0; i < n; i++)  Cboot.push_back(D[samples[i]]);

      cout << "Computing " << i << "th bootstrapped estimator..." << endl;
      compute_bootstrapped_estimator(samples,cloud,i,dist,Cboot,Nsubs,gain);
      cout << "Done." << endl;

      char namei[100]; sprintf(namei,"%s_mapper_%d.txt",cloud,i); char dgi[100]; sprintf(dgi,"%s_ExDg_%d",cloud,i);
      char command[100];
      sprintf(command,"cat %s | ../persistence/filtgraph 1 | ../persistence/pers 1 >> %s",namei,dgi); system(command);
      sprintf(command,"../persistence/bottleneck_dist %s %s >> %s",dg,dgi,list); system(command);
      sprintf(command, "rm %s %s",namei, dgi); system(command);

    }

    delete [] dist;

    return 0;
}
