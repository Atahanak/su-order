#ifndef _STATS
#define _STATS

template<typename vid_t, typename C, typename T>
void get_prefixsum_degrees(vid_t * xadj, C n, T * prefix_sum){
  memset(prefix_sum, 0, sizeof(vid_t) * n);
  for(vid_t i = 0; i < n; i++){
    prefix_sum[xadj[i+1] - xadj[i]]++;
  }
  for(vid_t i = 0; i < n; i++){
    prefix_sum[i+1] += prefix_sum[i];
  }
  //for(vid_t i = 0; i < n+1; i++){
  //  std::cout << prefix_sum[i] << " ";
  //}
  //std::cout << std::endl;
}
#endif
