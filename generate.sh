#! /bin/bash
        #com-dblp_c,4
        #com-lj_lcc_c,4
        #com-orkut_lcc_c,4
        #soc-LiveJournal_lcc_c,5
GRAPHS=(
        hyperlink2012_lcc_c,7
        soc-sinaweibo_lcc_c,10
        twitter_rv_lcc_c,7
        com-friendster_lcc_c,11
       )

DATA_PATH=/home/amro/graphs/uedgelist/
EXECS=(
      metis_order.out 
      )
ALGOS=(
      metis-naive
      )
for G in ${GRAPHS[*]}
  do
  for EXEC in ${EXECS[*]}
    do
    for ALGO in ${ALGOS[*]}
      do
        IFS=',' read -ra ADDR <<< "$G"
        GG=${ADDR[0]}
        PARTS=${ADDR[1]}
        INPUT=${DATA_PATH}${GG}.graph
        OUTPUT=${DATA_PATH}${GG}_${ALGO}.graph
        echo "bin/${EXEC} ${INPUT} ${ALGO} ${OUTPUT} $PARTS"
        bin/${EXEC} ${INPUT} ${ALGO} ${OUTPUT} $PARTS
      done
    done
  done

