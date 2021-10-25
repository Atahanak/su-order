#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cmath>

#include "../src/ordering.cpp"
#include "../src/utils/csr.cpp"

using namespace std;

int main(int argc, char * argv[]){
  if(argc != 2){
    cout << "Use: exec filename" << endl;
  } 
  
  ifstream infile(argv[1], ios::in);
  unsigned int * xadj; 
  unsigned int * adj; 
  unsigned int * tadj; 
  unsigned int * is; 
  unsigned int n, m;
  read_graph(infile, xadj, adj, tadj, is, n, m);
  cout << "# vertices: " << n << " # edges: " << m << endl;
  unsigned int * sorted; 
  degree(xadj, adj, n, m, sorted);
  cout << "sorted[n] deg = " << xadj[sorted[n-1]+1] - xadj[sorted[n-1]] << endl; 
  return 0;
}
