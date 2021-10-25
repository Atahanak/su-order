

using namespace std;

bool sortedge(const pair<int,int> &a,
    const pair<int,int> &b) {
  if(a.first == b.first) {
    return (a.second < b.second);
  } else {
    return (a.first < b.first);
  }
}

template<typename vid_t, typename eid_t>
int read_graph(ifstream & infile, vid_t *& xadj, vid_t *& adj, vid_t *& tadj, vid_t *& is, vid_t & n, eid_t& m);

template<typename vid_t, typename eid_t>
int read_graph(ifstream & infile, vid_t *& xadj, vid_t *& adj, vid_t *& tadj, vid_t *& is,  vid_t & n, eid_t & m){
    if(infile.is_open()) {
      vid_t u, v;
      eid_t edges_read = 0;
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
      m = edges.size();
      cout << "No edges is " << m << endl;

      sort(edges.begin(), edges.end(), sortedge);
      edges.erase( unique( edges.begin(), edges.end() ), edges.end() );

      //allocate the memory
      xadj = new vid_t[n + 1];
      adj = new vid_t[m];
      tadj = new vid_t[m];
      is = new vid_t[m];

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

      //ofstream outfile_bin(binary_name, ios::out | ios::binary);
      //if(outfile_bin.is_open()) {
      //  outfile_bin.write((char*)(&n), sizeof(int));
      //  outfile_bin.write((char*)(&m), sizeof(int));
      //  outfile_bin.write((char*)xadj, sizeof(int) * (n + 1));
      //  outfile_bin.write((char*)adj, sizeof(int) * m);
      //  outfile_bin.write((char*)tadj, sizeof(int) * m);
      //  outfile_bin.write((char*)is, sizeof(int) * m);
      //}
    } else {
      cout << "The file does not exist " << endl;
      return 1;
    }
    return 0;
}

