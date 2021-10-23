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

#include "utils/stats.cpp"

using namespace std;

bool sortedge(const pair<int,int> &a,
    const pair<int,int> &b) {
  if(a.first == b.first) {
    return (a.second < b.second);
  } else {
    return (a.first < b.first);
  }
}


//[Two improved algorithms for envelope and wavefront reduction]
//[ORDERING SYMMETRIC SPARSE MATRICES FOR SMALL PROFILE AND WAVEFRONT]

//b is blockcount
void testOrder(int *xadj, int* adj, int n, int b, int* order) {
  int* iperm = new int[n];
  for(int i = 0; i < n; i++) iperm[order[i]] = i;

  int max_bw = 0;
  double mean_bw = 0;

  int bsize = n / b;

  int** density = new int*[b];
  for(int i = 0; i < b; i++) {
    density[i] = new int[b];
    memset(density[i], 0, sizeof(int) * b); 
  }

  for(int i = 0; i < n; i++) {
    int u = iperm[i];
    int bu = u / bsize;
    if(bu == b) bu--;
    for(int ptr = xadj[i]; ptr < xadj[i+1]; ptr++) {
      int v = iperm[adj[ptr]];
      int bw = abs(u - v);
      max_bw = max(max_bw, bw);
      mean_bw += bw;

      int bv = v / bsize;
      if(bv == b) bv--;
      density[bu][bv]++;      
    }
  }
  mean_bw = (mean_bw + 0.0f) / xadj[n];

  cout << "BW stats -- Mean bw: " << mean_bw << " " << "Max bw: " << max_bw << endl;

  double para_mean_bw = 0;
  mean_bw = 0;
  int fblocks = 0;
  cout << "Printing blocks ----------------------------------------------" << endl;
  for(int i = 0; i < b; i++) {
    for(int j = 0; j < b; j++) {
      cout << std::setprecision(2) << density[i][j] / (xadj[n] + .0f) << "\t";
      if(density[i][j] > 0) {
        fblocks++;
      }
      int bw = abs(i - j);
      mean_bw += bw * density[i][j];
    }
    cout << endl;
  }
  cout << "---------------------------------------------------------------" << endl;
  cout << "Block BW stats -- No full blocks: " << fblocks << " Block BW: " << (mean_bw + 0.0f) / xadj[n] << endl;
  cout << "---------------------------------------------------------------" << endl;
}


#define TPB 128
void lastvertex(int* xadj, int* adj, int nov, int* distance, int source, int &last, int &dlast) {
  int jv, cv, noFrontier, noFrontier_prev, fdist = 0;
  int nextFsize, tEdgesRemain, ptr, eptr;
  int next = 1;
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

int peripheral(int* xadj, int* adj, int n, int start, int* distance, int* Q) {
  int r = start, r2;
  int rlevel = -1;
  int qlevel = 0;

  while(rlevel != qlevel) {
    //cout << "Finding peripheral: current dist = " << qlevel << endl;; 
    rlevel = qlevel;

#ifdef PARBFS
    lastvertex(xadj, adj, n, distance, r, r2, qlevel);
    r = r2;
#else

    for(int i = 0; i < n; i++) distance[i] = -1;
    int qrp = 0, qwp = 0;
    distance[r] = 0; Q[qwp++] = r;

    while(qrp < qwp) {
      int u = Q[qrp++];
      for(int ptr = xadj[u]; ptr < xadj[u+1]; ptr++) {
        int v = adj[ptr];
        if(distance[v] == -1) {
          distance[v] = distance[u] + 1;
          Q[qwp++] = v;
        }
      }
    }

    qlevel = 0;
    for(int i = 0; i < qrp; i++) {
      if(qlevel < distance[Q[i]]) {
        qlevel = distance[Q[i]];
        r = Q[i];	
      }
    }
#endif
  }    
  return r;
}

void rcm_sequential(int *xadj, int* adj, int n, int* Q, int* Qp, int* distance) {
  int* V = new int[n]; for(int i = 0; i < n; i++) V[i] = 0;
  priority_queue<pair<int, int> > PQ;
  int qrp = 0, qwp = 0;
  int reverse = n-1;

  for(int i = 0; i < n; i++) {
    if(V[i] == 0) {
      if(xadj[i] == xadj[i+1]) {
        Q[reverse--] = i;
        V[i] = 1;
        continue;
      }

      // cout << i << endl;
      int perv = peripheral(xadj, adj, n, i, distance, Qp);      
      V[perv] = 1; Q[qwp++] = perv;

      while(qrp < qwp) {
        int u = Q[qrp++];
        for(int ptr = xadj[u]; ptr < xadj[u+1]; ptr++) {
          int v = adj[ptr];
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
  for(int i = 0; i < n/2; i++) {
    int t = Q[i];
    Q[i] = Q[n - i - 1];
    Q[n - i - 1] = t;
  }
}

int main(int argc, char** argv) {
  if(argc != 2) {
    cout << "Use: exec filename (edge list problem) " << endl;
    return 1;
  }

  char binary_name[1024];
  sprintf(binary_name, "%s.met.bin", argv[1]);

  ifstream infile_bin(binary_name, ios::in | ios::binary);
  int *xadj, m, *adj, *is, n, *tadj;
  if(infile_bin.is_open()) {
    infile_bin.read((char*)(&n), sizeof(int));
    infile_bin.read((char*)(&m), sizeof(int));
    xadj = new int[n + 1];
    infile_bin.read((char*)xadj, sizeof(int) * (n + 1));

    adj = new int[m];
    infile_bin.read((char*)adj, sizeof(int) * m);

    tadj = new int[m];
    infile_bin.read((char*)tadj, sizeof(int) * m);

    is = new int[m];
    infile_bin.read((char*)is, sizeof(int) * m);
  } else {
    ifstream infile(argv[1]);
    if(infile.is_open()) {
      int u, v, edges_read = 0;
      n = 0;

      vector< std::pair<int, int> > edges;
      //vertices are 0-based                                                                                                                                                                                                                
      while (infile >> u >> v) {
        if(u != v) {
          edges.push_back(std::pair<int, int>(u, v));
          edges.push_back(std::pair<int, int>(v, u));

          n = max(n, u);
          n = max(n, v);

          edges_read++;
        }
      }
      n++;
      cout << "No vertices is " << n << endl;
      cout << "No read edges " << edges_read << endl;

      sort(edges.begin(), edges.end(), sortedge);
      edges.erase( unique( edges.begin(), edges.end() ), edges.end() );

      //allocate the memory                                                                                                                                                                                                                 
      xadj = new int[n + 1];

      m = edges.size();
      adj = new int[m];
      tadj = new int[m];
      is = new int[m];
      cout << "No edges is " << m << endl;

      //populate adj and xadj
      memset(xadj, 0, sizeof(int) * (n + 1));
      int mt = 0;
      for(std::pair<int, int>& e : edges) {
        xadj[e.first + 1]++;
        is[mt] = e.first;
        adj[mt++] = e.second;
      }

      for(int i = 1; i <= n; i++) {
        xadj[i] += xadj[i-1];
      }

      for(int i = 0; i < m; i++) {
        tadj[i] = xadj[adj[i]]++;
      }
      for(int i = n; i > 0; i--) {
        xadj[i] = xadj[i-1];
      }
      xadj[0] = 0;

      ofstream outfile_bin(binary_name, ios::out | ios::binary);
      if(outfile_bin.is_open()) {
        outfile_bin.write((char*)(&n), sizeof(int));
        outfile_bin.write((char*)(&m), sizeof(int));
        outfile_bin.write((char*)xadj, sizeof(int) * (n + 1));
        outfile_bin.write((char*)adj, sizeof(int) * m);
        outfile_bin.write((char*)tadj, sizeof(int) * m);
        outfile_bin.write((char*)is, sizeof(int) * m);
      }
    } else {
      cout << "The file does not exist " << endl;
      return 1;
    }
  }
  cout << "here " << endl;
  for(int i = 0; i < n; i++) {
    for(int j = xadj[i]; j < xadj[i + 1]; j++) {
      if(i != adj[tadj[j]]) {
        cout << "problem: " << i << " " << j << " " << adj[j] << " " << tadj[j] <<  endl;
      }
    }
  }
  cout << "there " << endl;
  int max_deg = 0, min_deg = n, deg;
  int* degs = new int[n];
  memset(degs, 0, sizeof(int) * n);
  int* ps_degs = new int[n];
  memset(ps_degs, 0, sizeof(int) * n);
  get_prefixsum_degrees(xadj, n, ps_degs);

  for(int u = 0; u < n; u++) {
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


  int nthreads = omp_get_max_threads();
  cout << "Running with " << nthreads << " threads \n";

  int b = 6;

  cout << "---------------------------------------------" << endl;
  cout << "Testing natural order" << endl;
  cout << "---------------------------------------------" << endl;
  int* order = new int[n];
  for(int i = 0; i < n; i++) order[i] = i;
  testOrder(xadj, adj, n, b, order);

  cout << endl;
  cout << "---------------------------------------------" << endl;
  cout << "Testing random order" << endl;
  cout << "---------------------------------------------" << endl;
  std::random_shuffle (order, order + n);
  testOrder(xadj, adj, n, b, order);

  cout << endl;
  cout << "---------------------------------------------" << endl;
  cout << "Testing RCM order" << endl;
  cout << "---------------------------------------------" << endl;
  int* Qb = new int[n];  
  int* distance = new int[n];
  double startt = omp_get_wtime();
  rcm_sequential(xadj, adj, n, order, distance, Qb);
  double endt = omp_get_wtime();
  cout << "It took " << endt - startt << " secs." << endl;
  testOrder(xadj, adj, n, b, order);
  return 0;
}



