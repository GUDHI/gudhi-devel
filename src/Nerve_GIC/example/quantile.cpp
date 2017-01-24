#include "mapper.h"

bool myComp(const pair<double,double> & P1, const pair<double,double> & P2){
  if(P1.first == P2.first)  return (P1.second <= P2.second);
  else  return P1.first < P2.first;
}

int main(int argc, char** argv){

  if(argc <= 3 || argc >= 5){cout << "./quant <file> <F/F^{-1}> <value>" << endl; return 0;}

  char* const file = argv[1];
  bool token = atoi(argv[2]);
  double value = atof(argv[3]);
  vector<double> V;
  vector<pair<double,double> > Q;

  double v; ifstream input(file); string line;
  while(getline(input,line)){
    stringstream stream(line); stream >> v; V.push_back(v);
  }

  int N = V.size();

  if(!token){
    int count = 0;
    for(int i = 0; i < N; i++)
      if(V[i] <= value)
        count++;
    cout << count*1.0/N << endl;
  }
  else{
    for(int i = 0; i < N; i++){
      int counti = 0;
      for(int j = 0; j < N; j++)
        if(V[j] <= V[i])
          counti++;
      Q.push_back(pair<double,double>(counti*1.0/N,V[i]));
    }
    sort(Q.begin(),Q.end(),myComp);
    int ind = 0;
    while(ind < N && Q[ind].first <= value)  ind++;
    if(ind==N)  cout << Q[ind-1].second << endl;
    else  cout << Q[ind].second << endl;
  }

  return 0;
}
