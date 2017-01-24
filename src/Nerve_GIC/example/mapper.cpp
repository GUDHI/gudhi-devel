#include "mapper.h"
#define SHIFT 1
#define DELTA 1
#define RESOLUTION 1

int main(int argc, char** argv){

  if(argc <= 10 || argc >= 12){cout << "./mapper <cloud_file> <function:/coordinate:/eccentricity:> <func_file/number/matrix_file> <graph:/parameter:/percentage:>" <<\
                                       " <graph_file/delta> <VNE> <cover:/uniform:> <cover_file/resolution> <gain> <mask>" << endl; return 0;}

  char* const cloud = argv[1];
  char* const funct = argv[2];
  bool normalized = atoi(argv[6]);
  char* const graph = argv[4];
  int mask;
  char* const covering = argv[7];

  Cover I; AdjacencyMatrix G; Cloud C;
  MapperElements M;

  cout << "Reading input cloud from file " << cloud << "..." << endl;
  C = read_cloud(cloud);
  cout << "  Done." << endl;

  double r,g; char namefunc[100];

  if (strcmp(funct, "function:") == 0){
    char* const func = argv[3];
    cout << "Reading input filter from file " << func << "..." << endl;
    read_function_from_file(func,C); sprintf(namefunc,"%s",func);
  }
  if (strcmp(funct, "coordinate:") == 0){
    int number = atoi(argv[3]);
    cout << "Using coordinate " << number << " as a filter..." << endl;
    read_coordinate(number,C); sprintf(namefunc,"Coordinate %d",number);
  }
  if (strcmp(funct, "eccentricity:") == 0){
    char* const matrix = argv[3];
    cout << "Computing eccentricity with distance matrix " << matrix << "..." << endl;
    compute_eccentricity(C,matrix); sprintf(namefunc,"eccentricity");
  }
  cout << "  Done." << endl;

  if (strcmp(graph, "graph:") == 0){
    char* const graph_name = argv[5];
    cout << "Reading neighborhood graph from file " << graph_name << "..." << endl;
    G = build_neighborhood_graph_from_file(C,graph_name);
  }
  if (strcmp(graph, "percentage:") == 0){
    double delta = atof(argv[5]);
    char name1[100]; sprintf(name1, "%s_dist", cloud); char* const name2 = name1;
    cout << "Computing neighborhood graph with delta percentage = " << delta << "..." << endl;
    G = build_neighborhood_graph_from_scratch_with_percentage(C,delta,name2,normalized);
  }
  if (strcmp(graph, "parameter:") == 0){
    double delta = atof(argv[5]);
    char name1[100]; sprintf(name1, "%s_dist", cloud); char* const name2 = name1;
    cout << "Computing neighborhood graph with delta parameter = " << delta << "..." << endl;
    G = build_neighborhood_graph_from_scratch_with_parameter(C,delta,name2,normalized);
  }
  cout << "  Done." << endl;

  cout << "Sorting cloud..." << endl;
  sort(C.begin(),C.end());
  cout << "  Done." << endl;

  if(strcmp(covering,"cover:") == 0){
    char* const cover = argv[8]; mask = atoi(argv[9]);
    cout << "Reading user-defined cover from file " << cover << "..." << endl;
    I = Cover(cover); I.sort_covering();
    assert (I.intervals[0].first <= C.begin()->func && I.intervals[I.intervals.size()-1].second >= (C.end()-1)->func);
  }
  else{
    double resolution = atof(argv[8]); r = resolution;
    double gain = atof(argv[9]); mask = atoi(argv[10]); g = gain;
    cout << "Computing uniform cover with resolution " << resolution << " and gain " << gain << "..." << endl;
    I = Cover(C.begin()->func, (C.end()-1)->func, resolution, gain, SHIFT, RESOLUTION);
  }
  I.proper_value();
  /*for (int i = 0; i < I.res; i++)
    cout << "  " << I.intervals[i].first << " " << I.intervals[i].second << " " << I.value[i] << endl;
  cout << "  Done." << endl;*/

  map<Point,vector<int> > clusters; vector<int> dum; dum.clear();
  for(int i = 0; i < G.size(); i++)  clusters.insert(pair<Point,vector<int> >(C[i],dum)); vector<double> col(G.size());
  cout << "Computing Mapper nodes..." << endl;
  M = MapperElts(G,I,clusters,col);
  cout << "  Done." << endl;

  cout << "Computing intersections..." << endl;
  map<int,double> colors; colors.clear(); map<int,int> num; num.clear(); SparseGraph MG;
  if(!DELTA)  MG = build_mapper_graph(M,colors,num);
  else  MG = build_mapperDelta_graph(M,G,clusters,colors,num);
  cout << "  Done." << endl;

  /*cout << "Computing Nerve..." << endl;
  SimplicialComplex Nerve = compute_Nerve(C,clusters);
  for(int i = 0; i < Nerve.size(); i++){
    for(int j = 0; j < Nerve[i].size(); j++)  cout << Nerve[i][j] << " ";
    cout << endl;
  }*/

  /*cout << "Computing GIC..." << endl;
  SimplicialComplex GIC = compute_GIC(G,clusters);
  for(int i = 0; i < GIC.size(); i++){
    for(int j = 0; j < GIC[i].size(); j++)  cout << GIC[i][j] << " ";
    cout << endl;
  }*/

  //cout << "Checking intersection crossings..." << endl;
  //vector<pair<int,int> > error_pairs = check_intersection_crossing(M,G);
  //cout << "Done." << endl;

  cout << "Writing outputs..." << endl;

  char mapp[11] = "mapper.dot";
  char mappg[11] = "mapper.txt";
  ofstream graphic(mapp); graphic << "graph Mapper {" << endl;
  ofstream graphicg(mappg);

  double maxv, minv; maxv = numeric_limits<double>::min(); minv = numeric_limits<double>::max();
  for (map<int,double>::iterator iit = colors.begin(); iit != colors.end(); iit++){
    if(iit->second > maxv){maxv = iit->second;}
    if(minv > iit->second){minv = iit->second;}
  }

  int k = 0; vector<int> nodes; nodes.clear();
  for (SparseGraph::iterator iit = MG.begin(); iit != MG.end(); iit++){
    if(num[iit->first] > mask){
      nodes.push_back(iit->first);
      graphic << iit->first << "[shape=circle fontcolor=black color=black label=\""
              << iit->first /*<< ":" << num[iit->first]*/ << "\" style=filled fillcolor=\""
              << (1-(maxv-colors[iit->first])/(maxv-minv))*0.6 << ", 1, 1\"]" << endl;
      k++;
    }
  }
  int ke = 0;
  for (SparseGraph::iterator it = MG.begin(); it != MG.end(); it++)
    for (int j = 0; j < it->second.size(); j++)
      if(num[it->first] > mask && num[it->second[j]] > mask){
        graphic << "  " << it->first << " -- " << it->second[j] << " [weight=15];" << endl; ke++;}
  graphic << "}"; graphic.close();

  graphicg << cloud << endl;
  graphicg << namefunc << endl;
  graphicg << r << " " << g << endl;
  graphicg << k << " " << ke << endl; int kk = 0;
  for (vector<int>::iterator iit = nodes.begin(); iit != nodes.end(); iit++){
    graphicg << kk << " " << colors[*iit] << " " << num[*iit] << endl; kk++;}
  for (SparseGraph::iterator it = MG.begin(); it != MG.end(); it++)
    for (int j = 0; j < it->second.size(); j++)
      if(num[it->first] > mask && num[it->second[j]] > mask)
        graphicg << it->first << " " << it->second[j] << endl;
  graphicg.close();

  cout << "  Done." << endl;

  char command[100];
  sprintf(command,"neato %s -Tpdf -o mapper.pdf",mapp); system(command);
  sprintf(command,"python visu.py"); system(command);
  sprintf(command,"firefox mapper_visualization_output.html"); system(command);
  sprintf(command,"evince mapper.pdf"); system(command);
  sprintf(command,"rm %s %s mapper_visualization_output.html mapper.pdf",mapp,mappg); system(command);

  /*
  if (error_pairs.size() > 0)
    for (vector<pair<int,int> >::iterator it = error_pairs.begin(); it != error_pairs.end(); it++)
      if( find(nodes.begin(),nodes.end(),it->first)!=nodes.end() && find(nodes.begin(),nodes.end(),it->second)!=nodes.end() )
        graphic << "  " << it->first << " -- " << it->second << " [weight=15];" << endl;
  */

  int dim = C[0].coord.size();
  if (dim <= 3){
    plotNG(dim,G);
    plotCover(dim,G,col,I);
  }

return 0;
}
