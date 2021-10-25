template <typename vid_t, typename eid_t>
void degree(vid_t * xadj, vid_t * adj, vid_t n, eid_t m, vid_t *& sorted);
template <typename vid_t, typename eid_t, typename w_t>
void rabbit(vid_t * xadj, eid_t * adj, w_t * w, vid_t * dendogram);


template <typename vid_t, typename eid_t>
void degree(vid_t * xadj, vid_t * adj, vid_t n, eid_t m, vid_t *& sorted){
  vid_t * counts = new vid_t[n];
  for(vid_t u = 0; u < n; u++){
    counts[xadj[u+1] - xadj[u]+1]++;
  }
  for(vid_t u = 1; u < n; u++){
    counts[u] += counts[u - 1];
  }
  sorted = new vid_t[n];
  memset(sorted, -1, sizeof(vid_t) * n);
  vid_t * mr = new vid_t[n]();
  for(vid_t u = 1; u < n; u++){
    vid_t ec = counts[xadj[u+1] - xadj[u]];
    sorted[ec + mr[ec]] = u;
    mr[ec]++;
  }
  delete mr;
  delete counts;
}

/* 
 * rabbit order -> https://www.rd.ntt/_assets/pdf/sic/team_researchers/other/araij2016ipdps.pdf 
 * Rabbit Order: Just-in-time Parallel Reordering for Fast Graph Analysis
 */
template <typename vid_t, typename eid_t, typename w_t>
void rabbit(vid_t * xadj, eid_t * adj, w_t * w, vid_t * dendogram){

}
