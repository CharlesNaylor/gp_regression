#!/bin/bash
for i in {1..5}
do
  ./gp_forecasts sample \
  num_samples=500 num_warmup=500 \
  random seed=8 \
  id=$i data file=../data/weekly/gp_$1.dat \
  output file=samples_$1_$i &
done

