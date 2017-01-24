#include "utils.h"

typedef map<Point, vector<Point> > AdjacencyMatrix; //sparse matrix of adjacencies
typedef map<int, vector<Point> > ConnectedComp; //a cc is characterized by its ID and the list of its points
typedef map<int,vector<int> > SparseGraph;
typedef map<int, pair<double,ConnectedComp> > MapperElements;
typedef vector<vector<int> > SimplicialComplex;

bool comp(const pair<double, double> & I1, const pair<double, double> & I2){return I1.second < I2.second;}

bool compSC(const vector<int> & V1, const vector<int> & V2){
  if (V1.size() != V2.size())  return V1.size() < V2.size();
  else{
    int tok = 0;
    while(V1[tok] == V2[tok] && tok < V1.size())  tok++;
    if(tok != V1.size())  return (V1[tok] < V2[tok]);
    else return 0;
  }
}

class Cover{

  public:

    int res;
    vector<pair<double, double> > intervals;
    vector<double> value;

    Cover(){}

    Cover(const double & minf, const double & maxf, const double & resolution, const double & gain, const double & shift, const bool & token){
      if(!token){
        double incr = (maxf-minf)/resolution; double x = minf; double alpha = (incr*gain)/(2-2*gain);
        double y = minf + shift*incr + alpha; pair<double, double> interm(x,y); intervals.push_back(interm);
        for(int i = 1; i < resolution-1; i++){
          x = minf + i*incr - alpha;
          y = minf + (i+1)*incr + alpha;
          pair<double, double> inter(x,y); intervals.push_back(inter);
        }
        x = minf + (resolution-1)*incr - alpha; y = maxf;
        pair<double, double> interM(x,y); intervals.push_back(interM); res = intervals.size();
      }
      else{
        double x = minf; double y = x + shift*resolution;
        while(y <= maxf && maxf - (y-gain*resolution) >= resolution){
          pair<double, double> inter(x,y); intervals.push_back(inter);
          x = y - gain*resolution;
          y = x + resolution;
        }
        pair<double, double> interM(x,maxf); intervals.push_back(interM); res = intervals.size();
      }
    }

    Cover(char* const & name){
      ifstream input(name); intervals.clear();
      if(input){
        string line;
        while(getline(input,line)){
          pair<double, double> inter; double x = numeric_limits<double>::max();
          stringstream stream(line);
          stream >> x; inter.first = x; stream >> x; inter.second = x;
          if(x != numeric_limits<double>::max()){assert (inter.first <= inter.second); intervals.push_back(inter);}
        }
        res = intervals.size();
      }
      else{cout << "  Failed to read file " << name << endl;}
    }

    void sort_covering(){sort(intervals.begin(),intervals.end(), comp);}

    void proper_value(){
      value.push_back(0.5*(intervals[0].first+intervals[1].first));
      for(int i = 1; i < this->res-1; i++)  value.push_back(0.5*(intervals[i-1].second+intervals[i+1].first));
      value.push_back(0.5*(intervals[res-2].second+intervals[res-1].second));
    }
};

bool check_inter(vector<Point> & v1, vector<Point> & v2){
  vector<Point> v3(v1.size()+v2.size());
  set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),v3.begin());
  return (v3[0].ID != -1 );
}




void compute_eccentricity(Cloud & C, char* const & name){

  int num_pts = C.size();
  double* ecc = new double[num_pts]; for(int i = 0; i < num_pts; i++)  ecc[i] = 0;
  ifstream input(name, ios::binary);
  if(input.good()){
    cout << "  Reading distances..." << endl;
    double d;
    for(int i = 0; i < num_pts; i++){
      for (int j = 0; j < num_pts; j++){
        //input >> d;
        input.read((char*) &d,8);
        if(d >= ecc[i])  ecc[i] = d;
      }
    }
    input.close();
  }
  else{cout << "No distance matrix found" << endl; return;}

  for(int i = 0; i < num_pts; i++)  C[i].func = ecc[i];
  delete [] ecc;
  return;
}

pair<double,double> compute_delta_from_file(const Cloud & C, char* const & name){

  ifstream input(name); double maximum = 0; double lip = 0;
  if(input){
    string line; int k = 0;
    while(getline(input,line)){
      int v = numeric_limits<int>::max(); stringstream stream(line);
      while(stream >> v){double dis = C[k].EuclideanDistance(C[v]); maximum = max(maximum, dis); lip = max(lip,abs(C[k].func-C[v].func)/dis);}
      k++;
    }
    return pair<double,double>(maximum,lip);
  }
  else{cout << "Failed to read file " << name << endl; return pair<double,double>(0,0);}

}

AdjacencyMatrix build_neighborhood_graph_from_file(const Cloud & C, char* const & name){

  int nb = C.size(); AdjacencyMatrix G; G.clear(); vector<Point> adj; adj.clear();
  for(int i = 0; i < nb; i++)
    G.insert(pair<Point, vector<Point> >(C[i],adj));
  ifstream input(name);
  if(input){
    string line; int num_edges = 0; int k = 0;
    while(getline(input,line)){
      vector<Point> ad; ad.clear(); int v = numeric_limits<int>::max(); stringstream stream(line);
      while(stream >> v)  ad.push_back(C[v]);
      if(v != numeric_limits<int>::max()){G[C[k]] = ad; k++; num_edges += ad.size();}
    }
    cout << "  " << 100*((double) num_edges/2)/((double) nb*(nb-1)/2) << "% of pairs selected." << endl;
    return G;
  }
  else{cout << "Failed to read file " << name << endl; return G;}

}

AdjacencyMatrix build_neighborhood_graph_from_scratch_with_percentage(const Cloud & C, const double & delta, char* const & name, const bool & VNE){

  ifstream input(name, ios::out | ios::binary);

  if(input.good()){
    cout << "  Reading distances..." << endl;
    int nb = C.size(); AdjacencyMatrix G; G.clear(); vector<Point> adj;
    double d; vector<vector<double> > dist; dist.clear();

    double m = 0;
    for(int i = 0; i < nb; i++){
      vector<double> dis; dis.clear();
      for (int j = 0; j < nb; j++){
        input.read((char*) &d,8); dis.push_back(d);
        if(m <= d)  m = d;
      }
      dist.push_back(dis);
    }
    input.close();

    cout << "  Done." << endl << "  Computing neighborhood graph...  ";

    cout.flush(); int num_edges = 0;
    for(int i = 0; i < nb; i++){
      adj.clear();
      for (int j = 0; j < nb; j++)
        if(dist[i][j] <= delta*m){adj.push_back(C[j]); num_edges++;}
      G.insert(pair<Point, vector<Point> >(C[i],adj));
    }
    cout << 100*((double) num_edges/2)/((double) nb*(nb-1)/2) << "% of pairs selected." << endl;
    return G;
  }
  else{
      cout << "  Computing distances..." << endl; vector<double> V; double m = 0;
      input.close(); ofstream output(name, ios::out | ios::binary);
      int nb = C.size(); AdjacencyMatrix G; G.clear(); vector<Point> adj;
      double d; vector<vector<double> > dist; dist.clear();

      if(VNE){
        int dim = C[0].coord.size(); int N = C.size();
        for(int i = 0; i < dim; i++){
          double meani = 0;
          for (int j = 0; j < nb; j++)  meani += C[j].coord[i]/N;
          double vari = 0;
          for (int j = 0; j < nb; j++)  vari += pow(C[j].coord[i]-meani,2)/N;
          V.push_back(vari);
        }
      }

      for(int i = 0; i < nb; i++){
        if( (int) floor( 100*((double) i)/((double) nb)+1 ) %10 == 0  ){cout << "\r" << floor( 100*((double) i)/((double) nb) +1) << "%" << flush;}
        vector<double> dis; dis.clear();
        for (int j = 0; j < nb; j++){
          if(!VNE){d = C[i].EuclideanDistance(C[j]); dis.push_back(d); if(m <= d)  m = d;}
          else{d = C[i].VarianceNormalizedEuclideanDistance(C[j],V); dis.push_back(d); if(m <= d)  m = d;}
          output.write((char*) &d,8);
        }
        dist.push_back(dis);
      }
      output.close();

    cout << endl << "  Done." << endl << "  Computing neighborhood graph...  ";

    cout.flush(); int num_edges = 0;
    for(int i = 0; i < nb; i++){
      adj.clear();
      for (int j = 0; j < nb; j++)
        if(dist[i][j] <= delta*m){adj.push_back(C[j]); num_edges++;}
      G.insert(pair<Point, vector<Point> >(C[i],adj));
    }
    cout << 100*((double) num_edges)/pow(nb,2) << "% of pairs selected." << endl;
    return G;
  }
}



AdjacencyMatrix build_neighborhood_graph_from_scratch_with_parameter(const Cloud & C, const double & delta, char* const & name, const bool & VNE){

    ifstream input(name, ios::out | ios::binary);

    if(input.good()){
      cout << "  Reading distances..." << endl;
      int nb = C.size(); AdjacencyMatrix G; G.clear(); vector<Point> adj;
      double d; vector<vector<double> > dist; dist.clear();

      for(int i = 0; i < nb; i++){
        vector<double> dis; dis.clear();
        for (int j = 0; j < nb; j++){
          input.read((char*) &d,8); dis.push_back(d);
        }
        dist.push_back(dis);
      }
      input.close();

      cout << "  Done." << endl << "  Computing neighborhood graph...  ";

      cout.flush(); int num_edges = 0;
      for(int i = 0; i < nb; i++){
        adj.clear();
        for (int j = 0; j < nb; j++)
          if(dist[i][j] <= delta && j != i){adj.push_back(C[j]); num_edges++;}
        G.insert(pair<Point, vector<Point> >(C[i],adj));
      }
      cout << 100*((double) num_edges/2)/((double) nb*(nb-1)/2) << "% of pairs selected." << endl;
      return G;
    }
    else{
      cout << "  Computing distances..." << endl; vector<double> V;
      input.close(); ofstream output(name, ios::out | ios::binary);
      int nb = C.size(); AdjacencyMatrix G; G.clear(); vector<Point> adj;
      double d; vector<vector<double> > dist; dist.clear();

      if(VNE){
        int dim = C[0].coord.size(); int N = C.size();
        for(int i = 0; i < dim; i++){
          double meani = 0;
          for (int j = 0; j < nb; j++)  meani += C[j].coord[i]/N;
          double vari = 0;
          for (int j = 0; j < nb; j++)  vari += pow(C[j].coord[i]-meani,2)/N;
          V.push_back(vari);
        }
      }

      for(int i = 0; i < nb; i++){
        if( (int) floor( 100*((double) i)/((double) nb)+1 ) %10 == 0  ){cout << "\r" << floor( 100*((double) i)/((double) nb) +1) << "%" << flush;}
        vector<double> dis; dis.clear();
        for (int j = 0; j < nb; j++){
          if(!VNE){d = C[i].EuclideanDistance(C[j]); dis.push_back(d);}
          else{d = C[i].VarianceNormalizedEuclideanDistance(C[j],V); dis.push_back(d);}
          output.write((char*) &d,8);
        }
        dist.push_back(dis);
      }
      output.close();

      cout << endl << "  Done." << endl << "  Computing neighborhood graph...  ";

      cout.flush(); int num_edges = 0;
      for(int i = 0; i < nb; i++){
        adj.clear();
        for (int j = 0; j < nb; j++)
          if(dist[i][j] <= delta && j != i){adj.push_back(C[j]); num_edges++;}
        G.insert(pair<Point, vector<Point> >(C[i],adj));
      }
      cout << 100*((double) num_edges)/pow(nb,2) << "% of pairs selected." << endl;
      return G;
    }
}




void dfs(AdjacencyMatrix & G, const Point & p, vector<Point> & cc, map<Point,bool> & visit){
  cc.push_back(p);
  visit[p] = true; int neighb = G[p].size();
  for (int i = 0; i < neighb; i++)
    if (  visit.find(G[p][i]) != visit.end() )
      if(  !(visit[G[p][i]])  )
        dfs(G,G[p][i],cc,visit);
}

ConnectedComp count_cc(AdjacencyMatrix & G, int & id, map<Point,vector<int> > & clusters){
  map<Point,bool> visit;
  for(AdjacencyMatrix::iterator it = G.begin(); it != G.end(); it++)
    visit.insert(pair<Point,bool>(it->first, false));
  ConnectedComp CC; CC.clear();
  if (!(G.empty())){
    for(AdjacencyMatrix::iterator it = G.begin(); it != G.end(); it++){
      if (  !(visit[it->first])  ){
        vector<Point> cc; cc.clear();
        dfs(G,it->first,cc,visit); int cci = cc.size();
        for(int i = 0; i < cci; i++)  clusters[cc[i]].push_back(id);
        CC.insert(pair<int, vector<Point> >(id++,cc));
      }
    }
  }
  return CC;
}

MapperElements MapperElts(AdjacencyMatrix & G, const Cover & I, map<Point,vector<int> > & clusters, vector<double> & col){

  map<int, pair<double,ConnectedComp> > mapper_elts;
  AdjacencyMatrix::iterator pos = G.begin();
  int id = 0;

  for(int i = 0; i < I.res; i++){

    AdjacencyMatrix prop; prop.clear();
    pair<double, double> inter1 = I.intervals[i];
    AdjacencyMatrix::iterator tmp = pos;

    if(i != I.res-1){
      if(i != 0){
        pair<double, double> inter3 = I.intervals[i-1];
        while(tmp->first.func < inter3.second && tmp != G.end()){
          col[tmp->first.ID] = 0;
          prop.insert(*tmp);
          tmp++;
        }
      }
      pair<double, double> inter2 = I.intervals[i+1];
      while(tmp->first.func < inter2.first && tmp != G.end()){
        col[tmp->first.ID] = i+1;
        prop.insert(*tmp);
        tmp++;
      }
      pos = tmp;
      while(tmp->first.func < inter1.second && tmp != G.end()){
        prop.insert(*tmp);
        tmp++;
      }

    }
    else{
      pair<double, double> inter3 = I.intervals[i-1];
      while(tmp->first.func < inter3.second && tmp != G.end()){
        col[tmp->first.ID] = 0;
        prop.insert(*tmp);
        tmp++;
      }
      while(tmp != G.end()){
        col[tmp->first.ID] = i+1;
        prop.insert(*tmp);
        tmp++;
      }

    }
    ConnectedComp CC = count_cc(prop,id,clusters);
    mapper_elts.insert(pair<int, pair<double,ConnectedComp> >(i, pair<double,ConnectedComp>(I.value[i],CC)));

  }

  return mapper_elts;
}

void find_all_simplices(vector<vector<int> > & lc, SimplicialComplex & SC, int & tok, vector<int> & simpl){
  int num_nodes = lc.size();
  if(tok == num_nodes-1){
    int num_clus = lc[tok].size();
    for(int i = 0; i < num_clus; i++){
      vector<int> simplex = simpl; simplex.push_back(lc[tok][i]);
      sort(simplex.begin(), simplex.end()); vector<int>::iterator iter = unique(simplex.begin(),simplex.end()); simplex.resize(distance(simplex.begin(),iter));
      SC.push_back(simplex); int dimension = simplex.size();
      if(dimension >= 2){for(int j = 0; j < dimension; j++){vector<int> face = simplex; face.erase(face.begin()+j); SC.push_back(face);}}
    }
  }
  else{
    int num_clus = lc[tok].size();
    for(int i = 0; i < num_clus; i++){
      vector<int> simplex = simpl; simplex.push_back(lc[tok][i]);
      find_all_simplices(lc, SC, ++tok, simplex);
    }
  }
}

void bron_kerbosch_II(vector<Point> R, vector<Point> P, vector<Point> X, AdjacencyMatrix & G, vector<vector<Point> > & cliques){

  if(P.size() == 0 && X.size() == 0){cliques.push_back(R); /*cout << "clique "; for(int i = 0; i < R.size(); i++)  cout << R[i].ID+1 << "  "; cout << endl;*/}
  else{
    if(P.size() != 0){
      Point pivot = P[0]; vector<Point> neighbors_pivot = G[pivot];

      sort(P.begin(),P.end()); //cout << "P = "; for(int i = 0; i < P.size(); i++)  cout << P[i].ID+1 << "  "; cout << endl;
      sort(neighbors_pivot.begin(),neighbors_pivot.end()); //cout << "pivot = " << pivot.ID+1 << endl;

      vector<Point> list(neighbors_pivot.size()+P.size()); vector<Point>::iterator it;
      it = set_difference(P.begin(),P.end(),neighbors_pivot.begin(),neighbors_pivot.end(),list.begin());
      list.resize(distance(list.begin(),it)); int nb = list.size();
      //cout << "N(pivot) = "; for(int i = 0; i < neighbors_pivot.size(); i++)  cout << G[pivot][i].ID+1 << "  "; cout << endl;
      //cout << "P minus N(pivot) = "; for(int i = 0; i < list.size(); i++)  cout << list[i].ID+1 << "  "; cout << endl;

      for(int i = 0; i < nb; i++){

        Point v = list[i]; //cout << "v = " << v.ID+1 << endl;
        vector<Point> neighbors_v = G[v]; sort(neighbors_v.begin(),neighbors_v.end());
        vector<Point> newR = R; newR.push_back(v);

        //cout << "P = "; for(int j = 0; j < P.size(); j++)  cout << P[j].ID+1 << "  "; cout << endl;
        //cout << "N(v) = "; for(int j = 0; j < neighbors_v.size(); j++)  cout << neighbors_v[j].ID+1 << "  "; cout << endl;

        vector<Point> newP(P.size()+neighbors_v.size()); it = set_intersection(P.begin(),P.end(),neighbors_v.begin(),neighbors_v.end(),newP.begin());
        newP.resize(it-newP.begin()); //cout << "new P = "; for(int j = 0; j < newP.size(); j++)  cout << newP[j].ID+1 << "  "; cout << endl;

        vector<Point> newX(X.size()+neighbors_v.size()); it = set_intersection(X.begin(),X.end(),neighbors_v.begin(),neighbors_v.end(),newX.begin());
        newX.resize(it-newX.begin()); //cout << "new X = "; for(int j = 0; j < newX.size(); j++)  cout << newX[j].ID+1 << "  "; cout << endl;

        bron_kerbosch_II(newR,newP,newX,G,cliques);
        it = find(P.begin(),P.end(),v); P.erase(it);
        X.push_back(v);

      }
    }
  }
}

SimplicialComplex compute_GIC(AdjacencyMatrix & G, map<Point,vector<int> > & clusters){
  SimplicialComplex SC; AdjacencyMatrix H;
  vector<Point> R, P, X; R.clear(); X.clear();
  for(AdjacencyMatrix::iterator it = G.begin(); it!=G.end(); it++){
    int numneigh = it->second.size(); int numclus = clusters[it->first].size();
    if(numclus > 1){H.insert(*it); P.push_back(it->first);}
    else{
      int currentcluster = clusters[it->first][0]; vector<Point> vectneigh;
      for(int i = 0; i < numneigh; i++){
        if(clusters[it->second[i]].size() > 1 || clusters[it->second[i]][0] != currentcluster)
          vectneigh.push_back(it->second[i]);
      }
      if(vectneigh.size() > 0){
        H.insert(pair<Point,vector<Point> >(it->first, vectneigh));
        P.push_back(it->first);
      }
    }
    vector<int> vertex(1); for(int i = 0; i < numclus; i++){vertex[0] = clusters[it->first][i]; SC.push_back(vertex);}
  }
  vector<vector<Point> > cliques; cout << "  Computing maximal cliques..." << endl;
  bron_kerbosch_II(R,P,X,H,cliques); cout << "  Done" << endl;
  int num_cliques = cliques.size();
  for(int i = 0; i < num_cliques; i++){
    int num_nodes = cliques[i].size(); vector<vector<int> > list_clusters(num_nodes);
    for(int j = 0; j < num_nodes; j++)  list_clusters[j] = clusters[cliques[i][j]];
    vector<int> simpl; int start = 0; find_all_simplices(list_clusters,SC,start,simpl);
  }
  sort(SC.begin(), SC.end(), compSC); SimplicialComplex::iterator iter = unique(SC.begin(), SC.end()); SC.resize(iter-SC.begin());
  return SC;
}

void add_simplices(vector<int> & clus, SimplicialComplex & SC){
  int dimension = clus.size();
  if(dimension == 1)  SC.push_back(clus);
  else{
    SC.push_back(clus);
    for(int i = 0; i < dimension; i++){
      vector<int> cl = clus; cl.erase(cl.begin()+i);
      add_simplices(cl,SC);
    }
  }
}

SimplicialComplex compute_Nerve(vector<Point> & C, map<Point,vector<int> > & clusters){
  SimplicialComplex SC; int num_pts = C.size();
  for(int i = 0; i < num_pts; i++){
    vector<int> clus = clusters[C[i]];
    add_simplices(clus,SC);
  }
  sort(SC.begin(), SC.end(), compSC); SimplicialComplex::iterator iter = unique(SC.begin(), SC.end()); SC.resize(iter-SC.begin());
  return SC;
}

SparseGraph build_mapper_graph(MapperElements & M, map<int, double> & colors, map<int,int> & num){

  SparseGraph MG; MG.clear(); vector<int> adj;
  int nb_elts = M.size(); adj.clear();

  #pragma omp parallel for
  for (int i = 0; i < nb_elts; i++){
    double c = M[i].first;
    for (ConnectedComp::iterator it = M[i].second.begin(); it != M[i].second.end(); it++){
      int p = it->first;
      #pragma omp critical
      {
        MG.insert(pair<int, vector<int> >(p,adj));
        int inside_nb = it->second.size();
        colors.insert(pair<int,double>(p,c)); num.insert(pair<int,int>(p,inside_nb));
      }
    }
  }


  for(int i = 0; i < nb_elts-1; i++){

    if(i==0){
      for (ConnectedComp::iterator it = M[i].second.begin(); it != M[i].second.end(); it++)
        sort(it->second.begin(),it->second.end());
    }

    for (ConnectedComp::iterator iit = M[i+1].second.begin(); iit != M[i+1].second.end(); iit++)
      sort(iit->second.begin(),iit->second.end());

    for (ConnectedComp::iterator it = M[i].second.begin(); it != M[i].second.end(); it++)
      for (ConnectedComp::iterator iit = M[i+1].second.begin(); iit != M[i+1].second.end(); iit++)
        if(check_inter(it->second, iit->second))
          MG[it->first].push_back(iit->first);
  }

  return MG;
}

SparseGraph build_mapperDelta_graph(MapperElements & M, AdjacencyMatrix & G, map<Point,vector<int> > & clusters, \
                                    map<int, double> & colors, map<int,int> & num){

  SparseGraph MG; MG.clear(); vector<int> adj;
  int nb_elts = M.size(); adj.clear();

  #pragma omp parallel for
  for (int i = 0; i < nb_elts; i++){
    double c = M[i].first;
    for (ConnectedComp::iterator it = M[i].second.begin(); it != M[i].second.end(); it++){
      int p = it->first;
      #pragma omp critical
      {
        MG.insert(pair<int, vector<int> >(p,adj));
        int inside_nb = it->second.size();
        colors.insert(pair<int,double>(p,c)); num.insert(pair<int,int>(p,inside_nb));
      }
    }
  }

  #pragma omp parallel for
  for(int i = 0; i < nb_elts-1; i++){
    double vali = M[i+1].first;
    for (ConnectedComp::iterator it = M[i].second.begin(); it != M[i].second.end(); it++){
      int sit = it->second.size(); int idit = it->first;
      for(int j = 0; j < sit; j++){
        Point p = it->second[j]; int neighp = G[p].size(); vector<int> clusp = clusters[p]; int numclus = clusp.size();
        for(int k = 0; k < numclus; k++)  if(clusp[k] > idit)  MG[idit].push_back(clusp[k]);
        for(int k = 0; k < neighp; k++){
          Point q = G[p][k]; int clusq = clusters[q].size();
          for(int l = 0; l < clusq; l++)
            if(clusters[q][l] > idit && colors[clusters[q][l]] == vali){
            #pragma omp critical
            {MG[idit].push_back(clusters[q][l]);}
            }
        }
      }
      sort(MG[idit].begin(),MG[idit].end()); vector<int>::iterator iter = unique(MG[idit].begin(),MG[idit].end());
      MG[idit].resize(distance(MG[idit].begin(),iter));
    }
  }

  return MG;
}

/*
vector<pair<int,int> > check_intersection_crossing(const vector<ConnectedComp> & M, AdjacencyMatrix & G){
  int nb_elts = M.size(); vector<pair<int,int> > pairs; pairs.clear();
  for (int i = 0; i < nb_elts-1; i+=2){
    ConnectedComp CC1 = M[i]; ConnectedComp CC2 = M[i+2];
    for(ConnectedComp::iterator it = CC1.begin(); it != CC1.end(); it++){
      int cc1_nb = it->second.size();
      for(ConnectedComp::iterator iit = CC2.begin(); iit != CC2.end(); iit++){
        int cc2_nb = iit->second.size(); int i = 0; bool point_in_inter = false; int count_interval_cross = 0;
        //cout << "  Checking nodes " << it->first << " and " << iit->first << "..." << endl;
        while(i < cc1_nb && !point_in_inter){
          Point p1 = it->second[i]; int j = 0;
          while(j < cc2_nb && !point_in_inter){
            Point p2 = iit->second[j];
            if(p1 == p2){point_in_inter = true;}
            if(  find(G[p1].begin(),G[p1].end(),p2)!=G[p1].end()  ){count_interval_cross++;}
            j++;
          }
          i++;
        }
        if(!point_in_inter && count_interval_cross > 0){
          cout << "  Warning: nodes " << it->first << " and " << iit->first << " separated whereas they should not: " << count_interval_cross << " intersection crossings" << endl;
          pairs.push_back(pair<int,int>(it->first,iit->first));
        }
      }
    }
  }
  return pairs;
}
*/

void plotNG(const int & dim, AdjacencyMatrix & G){
  char neigh[100]; int ned = 0; for (AdjacencyMatrix::iterator it = G.begin(); it != G.end(); it++){ned += it->second.size();}
  sprintf(neigh, "NG.vect");
  ofstream output(neigh);
  output << "VECT" << endl;
  output << ned + G.size() << " " << 2*ned + G.size() << " " << ned + G.size() << endl;
  for(int i = 0; i < ned; i++)  output << "2 ";
  for(int i = 0; i < G.size(); i++)  output << "1 ";
  output << endl;
  for(int i = 0; i < ned + G.size(); i++)  output << "1 ";
  output << endl;
  for (AdjacencyMatrix::iterator it = G.begin(); it != G.end(); it++){
    for(int j = 0; j < it->second.size(); j++){
      for(int k = 0; k < dim; k++)  output << it->first.coord[k] << " ";
      if(dim == 1)  output << "0 0  ";
      if(dim == 2)  output << "0  ";
      if(dim == 3)  output << "  ";
      for(int k = 0; k < dim; k++)  output << it->second[j].coord[k] << " ";
      if(dim == 1)  output << "0 0" << endl;
      if(dim == 2)  output << "0" << endl;
      if(dim == 3)  output << endl;
    }
  }
  for (AdjacencyMatrix::iterator it = G.begin(); it != G.end(); it++){
    for(int k = 0; k < dim; k++)  output << it->first.coord[k] << " ";
    if(dim == 1)  output << "0 0" << endl;
    if(dim == 2)  output << "0" << endl;
    if(dim == 3)  output << endl;
  }
  //for(int i = 0; i < ned + G.size(); i++){
  int plot_coord = 0;
  double c; double maxc = G.begin()->first.coord[2]; double minc = G.begin()->first.coord[plot_coord];
  if(dim == 3){
    for (AdjacencyMatrix::iterator it = G.begin(); it != G.end(); it++){
      if(it->first.coord[plot_coord] >= maxc)  maxc = it->first.coord[plot_coord];
      if(it->first.coord[plot_coord] <= minc)  minc = it->first.coord[plot_coord];
    }
  }
  for (AdjacencyMatrix::iterator it = G.begin(); it != G.end(); it++){
    if(dim == 3)  c = (it->first.coord[plot_coord]-minc)/(maxc-minc);
    else  c = 0;
    for(int j = 0; j < it->second.size(); j++){
      if(c<=0.5)  output << 1-2*c << " " << 2*c << " " << 0 << " 1" << endl;
      else  output << 0 << " " << -2*c+2 << " " << 2*c-1 << " 1" << endl;
    }
  }
  for (AdjacencyMatrix::iterator it = G.begin(); it != G.end(); it++){
    if(dim == 3)  c = (it->first.coord[plot_coord]-minc)/(maxc-minc);
    else  c = 0;
    if(c<=0.5)  output << 1-2*c << " " << 2*c << " " << 0 << " 1" << endl;
    else  output << 0 << " " << -2*c+2 << " " << 2*c-1 << " 1" << endl;
  }
  //   output << "0 0 1 1" << endl;
}

void plotCover(const int & dim, AdjacencyMatrix & G, const vector<double> & col, const Cover & I){
  int number = G.size();
  char neigh[100]; sprintf(neigh, "Cover.vect"); ofstream output(neigh);
  output << "VECT" << endl;
  output << number << " " << number << " " << number << endl;
  for(int i = 0; i < number; i++)  output << "1 ";
  output << endl;
  for(int i = 0; i < number; i++)  output << "1 ";
  output << endl;
  for (AdjacencyMatrix::iterator it = G.begin(); it != G.end(); it++){
    for(int k = 0; k < dim; k++)  output << it->first.coord[k] << " ";
    if(dim == 1)  output << "0 0" << endl;
    if(dim == 2)  output << "0" << endl;
    if(dim == 3)  output << endl;
  }
  int co = 1;
  for (AdjacencyMatrix::iterator it = G.begin(); it != G.end(); it++){
    double c = col[it->first.ID]/I.res;
    if(co != G.size()){
      if(c<=0.5){
        if(c==0)  output << 0 << " " << 0 << " " << 0 << " 1" << endl;
        else{
          if(c==1.0/I.res)  output << 1 << " " << 0 << " " << 0 << " 1" << endl;
          else  output << 1-2*c << " " << 2*c << " " << 0 << " 1" << endl;
        }
      }
      else  output << 0 << " " << -2*c+2 << " " << 2*c-1 << " 1" << endl;
      co++;}
    else{
      if(c<=0.5){
        if(c==0)  output << 0 << " " << 0 << " " << 0 << " 1" << endl;
        else{
          if(c==1.0/I.res)  output << 1 << " " << 0 << " " << 0 << " 1" << endl;
          else  output << 1-2*c << " " << 2*c << " " << 0 << " 1" << endl;
        }
      }
      else  output << 0 << " " << -2*c+2 << " " << 2*c-1 << " 1" << endl;
    }
  }
}
