#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <omp.h>
#include <random>
#include <cfloat>
#include <iomanip>
#include <queue>
#include <set>
#include <metis.h>
#ifdef __cplusplus
extern "C" {
#endif
#include <amd.h>
#ifdef __cplusplus
}
#endif


#include "utils/stats.cpp"

using namespace std;

template <typename T>
bool sortedge(const pair<T,T> &a,
    const pair<T,T> &b) {
  if(a.first == b.first) {
    return (a.second < b.second);
  } else {
    return (a.first < b.first);
  }
}

const string GRAPPOLO_BINARY_CONVERTER="grappolo/convertFileToBinary";
const string GRAPPOLO_CLUSTERING="grappolo/driverForGraphClustering";
enum MetisNumberingAlg {
  kNaive,
  kScoreCombined, // score each vertex by its total straddle edge count, and use score to greedily order vertices 
};

template <typename type>
void create_64bit_graph(unsigned int* xadj, unsigned int* adj, type*& new_xadj, type*& new_adj, long long n, long long m){
  new_xadj = new type[n+1];
  new_adj = new type[m];
  for (long long i = 0;i <n+1; i++) new_xadj[i] = xadj[i];
  for (long long i = 0;i <m; i++) new_adj[i] = adj[i];
}

struct count_name{
  unsigned long long count;
  unsigned int id;
};

struct count_name_compare{
  bool operator()(const count_name& lhs, const count_name& rhs){
    return  lhs.count < lhs.count;
  }
};

// Generates an order numbering using the following heuristic:
// For each part, the number of cut edges of each vertex is placed in an array
//  That array is made into a heap and elements are popped
//    For each vertex, if the number of cut edges to the left are bigger than the right, that vertex is added at the left, else its added at the right
template <typename T>
void score_metis_renumbering(T* part, unsigned int* order, long long n, int num_parts, unsigned int * xadj, unsigned int * adj){
  // count each part's size
  vector<long long> part_size(num_parts+1, 0);
  for (long long i = 0; i<n; i++) part_size[part[i]+1]++;
  long long max_part_size = 0;
  for (int i =0; i < num_parts;i++) max_part_size = max<long long>(max_part_size, part_size[i+1]);
  for (int i =0; i < num_parts;i++) part_size[i+1]+=part_size[i];
  count_name * boundary_counter = new count_name[max_part_size];
  long long* boundary_counter_sides = new long long[max_part_size*2];
  vector<vector<unsigned int>> names(num_parts);
  for (auto& x : names) x.reserve(max_part_size);
  for (long long i = 0; i<n; i++) names[part[i]].push_back(i);
  for (long long p= 0; p<num_parts; p++){
    memset(boundary_counter, 0, sizeof(count_name)*max_part_size);
    memset(boundary_counter_sides, 0, sizeof(long long)*max_part_size*2);
    for (int ui = 0; ui<names[p].size(); ui++){
      unsigned int u = names[p][ui];
      boundary_counter[ui].id = ui;
      for (int vi = xadj[u]; vi < xadj[u+1]; vi++){
        unsigned int partv = part[adj[vi]];
        if (partv != p){
          boundary_counter[ui].count++;
          if (partv < p) boundary_counter_sides[ui*2]++;
          else boundary_counter_sides[ui*2+1]++;
        }
      } 
    }
    make_heap(boundary_counter, boundary_counter+names[p].size(), count_name_compare());
    long long left = 0;
    long long right = names[p].size()-1;
    for (int ji = 0; ji < names[p].size(); ji++){
      auto boundary = boundary_counter[0];
      if (boundary_counter_sides[boundary.id*2] > boundary_counter_sides[boundary.id*2+1])
        order[names[p][boundary.id]] = part_size[p]+(left++);
      else
        order[names[p][boundary.id]] = part_size[p]+(right--);
      pop_heap(boundary_counter, boundary_counter+names[p].size()-ji, count_name_compare());
    }
  }
  delete[] boundary_counter;
  delete[] boundary_counter_sides;
}
template <typename T>
void naive_metis_renumbering(T* part, unsigned int* order, long long n, int num_parts){
  vector<long long> part_size(num_parts+1, 0);
  for (long long i = 0; i<n; i++) part_size[part[i]+1]++;
  for (int i =0; i < num_parts;i++) part_size[i+1]+=part_size[i];
  for (long long i =0; i<n; i++){
    order[i] = part_size[part[i]]++;
  }
}

bool metis_partitioning(unsigned int* xadj, unsigned int* adj, long long n, long long m, int num_parts, unsigned int * order, MetisNumberingAlg score_alg){
  idx_t* metis_xadj, *metis_adj;
  idx_t *nvtx = new idx_t(n);
  idx_t *ncon = new idx_t(1); // number of balancing constraints
  idx_t *nparts = new idx_t(num_parts);
  idx_t *objval = new idx_t;
  idx_t *part = new idx_t[*nvtx];
  idx_t options[METIS_NOPTIONS];
  double startt = omp_get_wtime();
  cout << "Creating 64 bit grap  ... " << flush;
  create_64bit_graph(xadj, adj, metis_xadj, metis_adj, n, m); 
  cout << "took " << omp_get_wtime()-startt << " secs." << endl;
  startt = omp_get_wtime();
  cout << "METIS partitioning  ... " << flush;
  METIS_PartGraphRecursive(nvtx, ncon, metis_xadj, metis_adj, NULL, NULL, NULL, nparts, NULL, NULL, NULL, objval, part);
  cout << "took " << omp_get_wtime()-startt << " secs." << endl;
  startt = omp_get_wtime();
  cout << "filling order array  ... " << flush;
  bool flag = true;
  if (score_alg == kNaive) 
    naive_metis_renumbering(part, order, n, num_parts);
  else if (score_alg == kScoreCombined){
    score_metis_renumbering(part, order, n, num_parts, xadj, adj);
  } else {
    cout << "Choose a legal renumbering algorithm for METIS\n";
    flag = false;
  }
  cout << "took " << omp_get_wtime()-startt << " secs." << endl;
  delete [] metis_xadj;
  delete [] metis_adj;
  delete nvtx;
  delete ncon; // number of balancing constraints
  delete nparts;
  delete objval;
  delete [] part;
  return flag;
}

bool grappolo_reordering(string filename, unsigned int*xadj, unsigned int* adj, long long n, unsigned int* order){
  
  // create the binary graph if it doesn't exist
  string binary_file_name = filename+".bin";
  ifstream bin_checker(binary_file_name);
  if (!bin_checker.good()){
    string command = GRAPPOLO_BINARY_CONVERTER +" -f 6 -o " +filename;
    char buffer[128];
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) cout << "popen failed\n";
    while (!feof(pipe)) {
      // use buffer to read and add to result
      if (fgets(buffer, 128, pipe) != NULL){
        cout << buffer;
      }
    }
    pclose(pipe);
  }
  // set the parameters required
  string command = GRAPPOLO_CLUSTERING+" -f 9 -o " + binary_file_name;
  // run the script
  long long num_clusters = -1;
  char buffer[128];
  FILE* pipe = popen(command.c_str(), "r");
  if (!pipe) cout << "popen failed\n";
  while (!feof(pipe)) {
    // use buffer to read and add to result
    if (fgets(buffer, 128, pipe) != NULL){
      cout << buffer << endl;
    }
    string l(buffer);
    if (l.find("Final number of clusters") != -1){
      stringstream ss(l);
      ss >> l >> l >> l >> l >> l >> num_clusters;
    }
  }
  if (num_clusters == -1){
    cout << "Couldn't find number of clusters\n";
    return false;
  }
  pclose(pipe);
  string cluster_info_file = binary_file_name+"_clustInfo";
  // read the clusterinfo file
  ifstream info(cluster_info_file);
  if (!info.good()){
    cout << "Failed to grappolo order\n";
    return false;
  }
  unsigned int* part = new unsigned int[n];
  char line[256];
  for (int i =0; i< n; i++){
    if (!info.getline(line, 256)){
      cout << "Not enough lines in cluster info file\n"; 
      delete [] part;
      return false;
    }
    stringstream is(line);
    is >> part[i];
  }
  
  // use one of the numbering functions to get numbers
  naive_metis_renumbering(part, order, n, num_clusters);

  return true;
} 



//[Two improved algorithms for envelope and wavefront reduction]
//[ORDERING SYMMETRIC SPARSE MATRICES FOR SMALL PROFILE AND WAVEFRONT]

//b is blockcount
void testOrder(unsigned int *xadj, unsigned int* adj, long long n, long long b, unsigned int* order) {
  //long long* iperm = new long long[n];
  //for(long long i = 0; i < n; i++) iperm[order[i]] = i;

  long long max_bw = 0;
  double mean_bw = 0;

  long long bsize = n / b;

  long long** density = new long long*[b];
  for(long long i = 0; i < b; i++) {
    density[i] = new long long[b];
    memset(density[i], 0, sizeof(long long) * b); 
  }

  for(long long i = 0; i < n; i++) {
    long long u = order[i];
    long long bu = u / bsize;
    if(bu == b) bu--;
    for(long long ptr = xadj[i]; ptr < xadj[i+1]; ptr++) {
      long long v = order[adj[ptr]];
      long long bw = abs(u - v);
      max_bw = max<long long>(max_bw, bw);
      mean_bw += bw;

      long long bv = v / bsize;
      if(bv == b) bv--;
      density[bu][bv]++;      
    }
  }
  mean_bw = (mean_bw + 0.0f) / xadj[n];

  cout << "BW stats -- Mean bw: " << mean_bw << " " << "Max bw: " << max_bw << endl;

  double para_mean_bw = 0;
  mean_bw = 0;
  long long fblocks = 0;
  cout << "Printing blocks \n----------------------------------------------" << endl;
  for(long long i = 0; i < b; i++) {
    for(long long j = 0; j < b; j++) {
      cout << std::setprecision(2) << density[i][j] / (xadj[n] + .0f) << "\t";
      if(density[i][j] > 0) {
        fblocks++;
      }
      long long bw = abs(i - j);
      mean_bw += bw * density[i][j];
    }
    cout << endl;
  }
  cout << "---------------------------------------------------------------" << endl;
  cout << "Block BW stats -- No full blocks: " << fblocks << " Block BW: " << (mean_bw + 0.0f) / xadj[n] << endl;
  cout << "---------------------------------------------------------------" << endl;
}


#define TPB 128
void lastvertex(unsigned int* xadj, unsigned int* adj, long long nov, long long* distance, long long source, long long &last, long long &dlast) {
  long long jv, cv, noFrontier, noFrontier_prev, fdist = 0;
  long long nextFsize, tEdgesRemain, ptr, eptr;
  long long next = 1;
  double alpha = 4, beta = 4;
  last = dlast = -1;
  next = 1;

#pragma omp parallel for
  for(jv = 0; jv < nov; jv++)  {
    distance[jv] = -1;
  }
  distance[source] = 0;

  noFrontier = 1;
  nextFsize = xadj[source + 1] - xadj[source];
  tEdgesRemain = xadj[nov] - nextFsize;

  while(noFrontier > 0) {
    noFrontier_prev = noFrontier;
    noFrontier = nextFsize = 0;

    if(next == 1) {
#pragma omp parallel for private(cv, ptr, eptr) reduction(+:nextFsize, noFrontier) schedule(dynamic, TPB)
      for(jv = 0; jv < nov; jv++)  {
        if(distance[jv] == fdist) {
          eptr = xadj[jv + 1];
          for(ptr = xadj[jv]; ptr < eptr; ptr++) {
            cv = adj[ptr];
            if(distance[cv] == -1) {
              noFrontier++;
              nextFsize += xadj[cv+1] - xadj[cv];
              distance[cv] = fdist + 1;
              last = cv;
              dlast = fdist + 1;
            }
          }
        }
      }
#ifdef TM
      printf("TD - %f\n", omp_get_wtime() - t1);
#endif
    } else if(next == 2) {
#pragma omp parallel for private(cv, ptr, eptr) reduction(+:nextFsize, noFrontier) schedule(dynamic, TPB)
      for(jv = 0; jv < nov; jv++)  {
        if(distance[jv] == -1) {
          eptr = xadj[jv + 1];
          for(ptr = xadj[jv]; ptr < eptr; ptr++) {
            cv = adj[ptr];
            if(distance[cv] == fdist) {
              distance[jv] = fdist + 1;
              last = jv;
              dlast = fdist + 1;
              noFrontier++;
              nextFsize += xadj[jv+1] - xadj[jv];
              break;
            }
          }
        }
      }
    }
    fdist++;
    tEdgesRemain -= nextFsize;

    //-------- transition logic -----------                                                                                                                                                                 
    if(next == 1) {
      if((nextFsize > (tEdgesRemain / alpha)) && (noFrontier > noFrontier_prev)) {
        next = 2;
      }
    } else if(next == 2) {
      if((noFrontier < (nov / beta)) && (noFrontier < noFrontier_prev)) {
        next = 1;
      }
    }
  }
}

long long peripheral(unsigned int* xadj, unsigned int* adj, long long n, long long start, long long* distance, long long* Q) {
  long long r = start, r2;
  long long rlevel = -1;
  long long qlevel = 0;

  while(rlevel != qlevel) {
    //cout << "Finding peripheral: current dist = " << qlevel << endl;; 
    rlevel = qlevel;

#ifdef PARBFS
    lastvertex(xadj, adj, n, distance, r, r2, qlevel);
    r = r2;
#else

    for(long long i = 0; i < n; i++) distance[i] = -1;
    long long qrp = 0, qwp = 0;
    distance[r] = 0; Q[qwp++] = r;

    while(qrp < qwp) {
      long long u = Q[qrp++];
      for(long long ptr = xadj[u]; ptr < xadj[u+1]; ptr++) {
        long long v = adj[ptr];
        if(distance[v] == -1) {
          distance[v] = distance[u] + 1;
          Q[qwp++] = v;
        }
      }
    }

    qlevel = 0;
    for(long long i = 0; i < qrp; i++) {
      if(qlevel < distance[Q[i]]) {
        qlevel = distance[Q[i]];
        r = Q[i];	
      }
    }
#endif
  }    
  return r;
}
template <typename T, typename L>
bool confirm_ordering(T* order, L n){
  set<T> seen;
  int duplicates = 0;
  for (int i =0; i < n; i++){
    if (seen.find(order[i]) != seen.end()){
      cout << "duplicate order " << i << " " << order[i] << "\n";
      duplicates++;
    }
    if (order[i]>=n){
      cout << "Order too big" << endl;
      return false;
    }
    seen.insert(order[i]);
  }
  if (seen.size() != n){
    cout << "Not every vertex has an order" << endl;
    return false;
  }
  if (duplicates>0){
    cout << "duplicates " << duplicates << endl;
    return false;
  }
  return true;
}
template <typename V, typename E, typename O, typename L>
bool confirm_renumbered_graph(V* xadj, V* renumbered_xadj, E* adj, E* renumbered_adj, O* order, O* inverse_order, L n){
  for (L i =0; i<n; i++){
    if (xadj[i+1]-xadj[i] != renumbered_xadj[order[i]+1]-renumbered_xadj[order[i]]){
      cout << "Mismatch number of edges\n";
      return false;
    }
    for (E edge = 0; edge < xadj[i+1]-xadj[i]; edge++){
      if (order[adj[xadj[i]+edge]]!=renumbered_adj[renumbered_xadj[order[i]]+edge]){
        cout << "Bad edge\n";
        return false;
      }
    }
  }
  return true;
}
template <typename V, typename E, typename C>
void print_graph(E* xadj, V* adj, C n, string filename){
  ofstream fout(filename);
  for (C i = 0; i < n; i++){
    sort(adj+xadj[i], adj+xadj[i+1]);
  }
  for (long long i =0; i<n; i++){
    for (long long j = xadj[i]; j < xadj[i+1]; j++){
      if (adj[j] < i) continue;
      fout << i << " " << adj[j] << '\n';
    }
  }
}
template <typename V, typename E, typename L, typename O>
void reorder_graph(E* xadj, V* adj, L n, O* order, E*new_xadj, V* new_adj, O* inverse_map){
  new_xadj[0] = 0;
  // an inverse map: iperm[i] = original id of vertex with new id i
  for (L i = 0; i < n; i++) inverse_map[order[i]] = i;
  for (L i = 0; i < n; i++){
    V old_name = inverse_map[i];
    L num_edges = xadj[old_name+1]-xadj[old_name];
    new_xadj[i+1] = new_xadj[i]+num_edges; 
    for (L j = 0; j < num_edges; j++){
      new_adj[new_xadj[i]+j] = order[adj[xadj[old_name]+j]];
    }
  }
}

bool amd_sequential(unsigned int *xadj, unsigned int* adj, long long n, unsigned int* order) {
  long *xadj_long, *adj_long;
  create_64bit_graph(xadj, adj, xadj_long, adj_long, n, xadj[n]);
  long * i_order = new long[n];
  double *Control = new double[AMD_CONTROL];
  double *Info = new double[AMD_INFO];
  int status = amd_l_order(n, xadj_long, adj_long, i_order, Control, Info);
  if (status != 0){
    cout << "AMD failed\n";
    return false;
  }
  for (int i =0; i< n; i++) order[i_order[i]]=i;
  return true;
}

void rcm_sequential(unsigned int *xadj, unsigned int* adj, long long n, unsigned int* Q, long long* Qp, long long* distance) {
  long long* V = new long long[n]; for(long long i = 0; i < n; i++) V[i] = 0;
  priority_queue<pair<long long, long long> > PQ;
  long long qrp = 0, qwp = 0;
  long long reverse = n-1;

  for(long long i = 0; i < n; i++) {
    if(V[i] == 0) {
      if(xadj[i] == xadj[i+1]) {
        Q[reverse--] = i;
        V[i] = 1;
        continue;
      }

      // cout << i << endl;
      long long perv = peripheral(xadj, adj, n, i, distance, Qp);      
      V[perv] = 1; Q[qwp++] = perv;

      while(qrp < qwp) {
        long long u = Q[qrp++];
        for(long long ptr = xadj[u]; ptr < xadj[u+1]; ptr++) {
          long long v = adj[ptr];
          if(V[v] == 0) {
            PQ.push(make_pair(xadj[v + 1] - xadj[v], v));
            V[v] = 1;
          }
        }

        while(!PQ.empty()) {
          Q[qwp++] = PQ.top().second;;
          PQ.pop();
        }
      }
    }
  }

  //Reverse
  for(long long i = 0; i < n/2; i++) {
    long long t = Q[i];
    Q[i] = Q[n - i - 1];
    Q[n - i - 1] = t;
  }
  for (long long i = 0; i < n; i++){
    Qp[i] = Q[i];
  }
  for (long long i = 0; i < n; i++){
    Q[Qp[i]] = i;
  }
}
vector<string> split(string input, string delimiter){
  size_t pos = 0;
  std::string token;
  vector<string>out;
  while ((pos = input.find(delimiter)) != std::string::npos) {
    token = input.substr(0, pos);
    out.push_back(token);
    input.erase(0, pos + delimiter.length());
  }
  out.push_back(input);
  return out;
}
const string METIS_NAME = "_METIS";
const string RCM_NAME = "_RCM";
vector<string> allowed_algorithms = {"rcm", "metis-naive", "metis-score", "random", "no-order", "grappolo", "amd"};

bool reorder_and_print_graph(unsigned int *xadj, unsigned int*adj, long long n, long long m, unsigned int* order, string filename){
  unsigned int* inverse_order = new unsigned int[n];
  for (int i =0; i<n; i++) inverse_order[order[i]] = i;
  unsigned int *new_xadj = new unsigned int[n+1];
  unsigned int *new_adj = new unsigned int[m];
  double startt = omp_get_wtime();
  cout << "Renumbering graph  ... " << flush;
  reorder_graph(xadj, adj, n, order, new_xadj, new_adj, inverse_order);
  cout << "took " << omp_get_wtime()-startt << " secs." << endl;
  startt = omp_get_wtime();
  cout << "Checking renumbering correctness  ... " << flush;
  bool confirmed = confirm_renumbered_graph(xadj, new_xadj, adj, new_adj, order, inverse_order, n);
  cout << "took " << omp_get_wtime()-startt << " secs." << endl;
  if (!confirmed){
    cout << "Graph renaming is wrong\n";
    delete [] inverse_order;
    delete [] new_xadj;
    delete [] new_adj;
    return false;
  }
  if (filename[0]!='-'){
    //auto graph_path_split = split(argv[1], "/");
    //auto graph_name_split = split(graph_path_split[graph_path_split.size()-1], ".graph");

    //string rcm_name = graph_name_split[0]+RCM_NAME+".graph";
    string name = filename;
    cout << "Printing the renamed graph to " << name << endl;
    startt = omp_get_wtime();
    cout << "Printing graph  ... " << flush;
    print_graph(new_xadj, new_adj, n, name);
    cout << "took " << omp_get_wtime()-startt << " secs." << endl;
  }
  delete [] inverse_order;
  delete [] new_xadj;
  delete [] new_adj;
  return true;

}
int main(int argc, char** argv) {
  if(argc < 4) {
    cout << "Use: exec filename algorithm output_filename [block_number] \n";
    cout << "(if you don't want to print out a reordered graph, pass - as the output_filename)\n ";
    cout << "Algorithms:\n:";
    for (auto alg : allowed_algorithms) cout << alg << endl;
    return 1;
  }

  char binary_name[1024];
  sprintf(binary_name, "%s.met.bin", argv[1]);
  string algorithm = argv[2];
  cout << "Using algorithm: " << algorithm << endl;
  if (std::find(allowed_algorithms.begin(), allowed_algorithms.end(), algorithm) == allowed_algorithms.end()){
    cout << "Please use one of the allowed algorithms:\n";
    for (auto alg : allowed_algorithms) cout << alg << endl;
    return 1;
  }
  
  // graph doesn't contain multi-edges and selfloops occur only once
  ifstream infile_bin(binary_name, ios::in | ios::binary);
  unsigned int *xadj, m, *adj, *is, n, *tadj;
  if(infile_bin.is_open()) {
    infile_bin.read((char*)(&n), sizeof(unsigned int));
    infile_bin.read((char*)(&m), sizeof(unsigned int));
    xadj = new unsigned int[n + 1];
    infile_bin.read((char*)xadj, sizeof(unsigned int) * (n + 1));

    adj = new unsigned int[m];
    infile_bin.read((char*)adj, sizeof(unsigned int) * m);

    tadj = new unsigned int[m];
    infile_bin.read((char*)tadj, sizeof(unsigned int) * m);

    is = new unsigned int[m];
    infile_bin.read((char*)is, sizeof(unsigned int) * m);
  } else {
    ifstream infile(argv[1]);
    if(infile.is_open()) {
      long long u, v, edges_read = 0;
      n = 0;

      vector< std::pair<long long, long long> > edges;
      //vertices are 0-based                                                                                                                                                                                                                
      while (infile >> u >> v) {
        if(u != v) {
          edges.push_back(std::pair<long long, long long>(u, v));
          edges.push_back(std::pair<long long, long long>(v, u));

          n = max<long long>(n, u);
          n = max<long long>(n, v);

          edges_read++;
        }
      }
      n++;
      cout << "No vertices is " << n << endl;
      cout << "No read edges " << edges_read << endl;

      sort(edges.begin(), edges.end(), sortedge<long long>);
      edges.erase( unique( edges.begin(), edges.end() ), edges.end() );

      //allocate the memory                                                                                                                                                                                                                 
      xadj = new unsigned int[n + 1];

      m = edges.size();
      adj = new unsigned int[m];
      tadj = new unsigned int[m];
      is = new unsigned int[m];
      cout << "No edges is " << m << endl;

      //populate adj and xadj
      memset(xadj, 0, sizeof(unsigned int) * (n + 1));
      long long mt = 0;
      for(std::pair<long long, long long>& e : edges) {
        xadj[e.first + 1]++;
        is[mt] = e.first;
        adj[mt++] = e.second;
      }

      for(long long i = 1; i <= n; i++) {
        xadj[i] += xadj[i-1];
      }

      for(long long i = 0; i < m; i++) {
        tadj[i] = xadj[adj[i]]++;
      }
      for(long long i = n; i > 0; i--) {
        xadj[i] = xadj[i-1];
      }
      xadj[0] = 0;

      ofstream outfile_bin(binary_name, ios::out | ios::binary);
      if(outfile_bin.is_open()) {
        outfile_bin.write((char*)(&n), sizeof(unsigned int));
        outfile_bin.write((char*)(&m), sizeof(unsigned int));
        outfile_bin.write((char*)xadj, sizeof(unsigned int) * (n + 1));
        outfile_bin.write((char*)adj, sizeof(unsigned int) * m);
        outfile_bin.write((char*)tadj, sizeof(unsigned int) * m);
        outfile_bin.write((char*)is, sizeof(unsigned int) * m);
      }
    } else {
      cout << "The file does not exist " << endl;
      return 1;
    }
  }
  for(long long i = 0; i < n; i++) {
    for(long long j = xadj[i]; j < xadj[i + 1]; j++) {
      if(i != adj[tadj[j]]) {
        cout << "problem: " << i << " " << j << " " << adj[j] << " " << tadj[j] <<  endl;
      }
    }
  }
  long long max_deg = 0, min_deg = n, deg;
  long long* degs = new long long[n];
  memset(degs, 0, sizeof(long long) * n);
  long long* ps_degs = new long long[n];
  memset(ps_degs, 0, sizeof(long long) * n);
  get_prefixsum_degrees(xadj, n, ps_degs);

  for(long long u = 0; u < n; u++) {
    deg = (xadj[u + 1] - xadj[u]);
    degs[deg]++;
    if(deg < min_deg) {min_deg = deg;}
    if(deg > max_deg) {max_deg = deg;}
  }

  cout << "---------------------------" << endl;
  cout << "No vertices is " << n << endl;
  cout << "No edges is " << m << endl;
  cout << "---------------------------" << endl;
  cout << "Min deg: " << min_deg << endl;
  cout << "Max deg: " << max_deg << endl;
  cout << "Avg deg: " << ((float)m)/n << endl;
  cout << "---------------------------" << endl;
  cout << "# deg 0: " << degs[0] << endl;
  cout << "# deg 1: " << degs[1] << endl;
  cout << "# deg 2: " << degs[2] << endl;
  cout << "# deg 3: " << degs[3] << endl;
  cout << "---------------------------" << endl;
  cout << "# deg>32: " << n-ps_degs[32] << endl;
  cout << "# deg>64: " << n-ps_degs[64] << endl;
  cout << "# deg>128: " << n-ps_degs[128] << endl;
  cout << "# deg>256: " << n-ps_degs[256] << endl;
  cout << "# deg>512: " << n-ps_degs[512] << endl;
  cout << "# deg>1024: " << n-ps_degs[1024] << endl;
  cout << "---------------------------" << endl << endl;


  long long nthreads = omp_get_max_threads();
  cout << "Running with " << nthreads << " threads \n";

  long long b = 6;
  if (argc > 4) b = atoi(argv[4]);

  unsigned int* order = new unsigned int[n];
  for(long long i = 0; i < n; i++) order[i] = i;
  /*
  cout << "---------------------------------------------" << endl;
  cout << "Testing natural order" << endl;
  cout << "---------------------------------------------" << endl;
  testOrder(xadj, adj, n, b, order);
  cout << endl;
  */
  if (algorithm == "random"){ 
    cout << "---------------------------------------------" << endl;
    cout << "Random order" << endl;
    cout << "---------------------------------------------" << endl;
    std::random_shuffle (order, order + n);
    confirm_ordering(order, n);
  }
  else if (algorithm == "rcm"){
    cout << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Testing RCM order" << endl;
    cout << "---------------------------------------------" << endl;
    long long* Qb = new long long[n];  
    long long* distance = new long long[n];
    double startt = omp_get_wtime();
    cout << "Generating RCM ordering ... " << flush;
    rcm_sequential(xadj, adj, n, order, distance, Qb);
    delete [] Qb;
    delete [] distance;
    double endt = omp_get_wtime();
    cout << "took " << endt - startt << " secs." << endl;
    confirm_ordering(order, n);
  } else if (algorithm == "metis-naive"){
    bool worked = metis_partitioning(xadj, adj, n, m, b, order, kNaive);

    if (!worked || !confirm_ordering(order, n)){
      cout << "Failed\n";
      return 1;
    }
  } else if (algorithm == "metis-score"){
    bool worked = metis_partitioning(xadj, adj, n, m, b, order, kScoreCombined);
    if (!worked || !confirm_ordering(order, n)){
      cout << "Failed\n";
      return 1;
    }
  } else if (algorithm == "grappolo"){
    grappolo_reordering(argv[1], xadj, adj, n, order);
  } else if (algorithm == "amd"){
    bool worked = amd_sequential(xadj, adj, n, order);
    if (!worked || !confirm_ordering(order, n)){
      cout << "Failed\n";
      return 1;
    }
  } else if (algorithm == "no-order"){}
  testOrder(xadj, adj, n, b, order);
  bool renumber_correct = reorder_and_print_graph(xadj, adj, n, m, order, argv[3]);
  if (!renumber_correct) {
    cout << "Renumbering the graph failed. Did not print\n";
  }
  
  return 0;
}



