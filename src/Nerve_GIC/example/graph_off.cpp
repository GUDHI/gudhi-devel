#include <iostream>
#include <set>
#include <fstream>
#include <vector>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <algorithm>

using namespace std;

int main (int argc, char **argv) {

  char* const fileoff = argv[1]; ifstream input(fileoff); char ngoff[100]; sprintf(ngoff,"%s_NG", fileoff); ofstream output(ngoff);
  char buf[256];
  vector<vector<int> > NG;

  input.getline(buf, 255);  // skip "OFF"
  int n, m;
  input >> n >> m;
  input.getline(buf, 255);  // skip "0"

  // read vertices
  double x,y,z; vector<int> ng; int nn = n;
  while(nn-->0) {
    input >> x >> z >> y; NG.push_back(ng);
  }

  // read triangles
  int d, p, q, s;
  while (m-->0) {
    input >> d >> p >> q >> s;
    NG[p].push_back(q); NG[p].push_back(s);
    NG[q].push_back(p); NG[q].push_back(s);
    NG[s].push_back(q); NG[s].push_back(p);
  }

  for(int i = 0; i < n; i++){
    vector<int> ng = NG[i];
    sort(ng.begin(),ng.end()); vector<int>::iterator iter = unique(ng.begin(),ng.end()); ng.resize(distance(ng.begin(),iter));
    int size = ng.size();
    for(int j = 0; j < size; j++)
      output << ng[j] << " ";
    output << endl;
  }

  return 0;

}
