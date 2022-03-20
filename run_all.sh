#!/bin/bash

DIR="results"
# change the path
DIR_data="/usr/ALENEX_datasets"
EXE_FILE="build/my_benchmark"

OPTIONS="--meva --plai --elia --rrrv --rlev --lave --lavo --encv --ds2i --s18v"

mkdir -p $DIR

time {
  $EXE_FILE --entr $DIR_data/DNA_1 $DIR_data/5GRAM_1 $DIR_data/URL_1 $DIR_data/GOV2/10M-/1_gov.bin.docs \
    >>$DIR/entropy_comparison.csv \
    2>>$DIR/comparison.err

  $EXE_FILE $OPTIONS $DIR_data/DNA_1 $DIR_data/DNA_2 $DIR_data/DNA_3 \
    $DIR_data/5GRAM_1 $DIR_data/5GRAM_2 $DIR_data/5GRAM_3 \
    $DIR_data/URL_1 $DIR_data/URL_2 $DIR_data/URL_3 \
    >>$DIR/DNA_5GRAM_URL_comparison.csv \
    2>>$DIR/comparison.err

  $EXE_FILE $OPTIONS $DIR_data/GOV2/100K-1M/* \
    >>$DIR/GOV2_100K-1M_comparison.csv \
    2>>$DIR/comparison.err

  $EXE_FILE $OPTIONS $DIR_data/GOV2/1M-10M/* \
    >>$DIR/GOV2_1M-10M_comparison.csv \
    2>>$DIR/comparison.err

  $EXE_FILE $OPTIONS $DIR_data/GOV2/10M-/* \
    >>$DIR/GOV2_10M-_comparison.csv \
    2>>$DIR/comparison.err

  NUM_COLUMNS=$(head -1 $DIR/GOV2_10M-_comparison.csv | sed 's/[^,]//g' | wc -c)

  printf "'GOV2_AVG_10M-'" >>$DIR/GOV2_averages.csv
  for i in $(seq 2 "$NUM_COLUMNS"); do
    printf "," >>$DIR/GOV2_averages.csv
    awk -v var="$i" -F ',' '{ total += $var } END { printf "%.3e", total/(NR-1) }' $DIR/GOV2_10M-_comparison.csv >>$DIR/GOV2_averages.csv
  done
  printf "\n" >>$DIR/GOV2_averages.csv

  printf "'GOV2_AVG_1M-10M'" >>$DIR/GOV2_averages.csv
  for i in $(seq 2 "$NUM_COLUMNS"); do
    printf "," >>$DIR/GOV2_averages.csv
    awk -v var="$i" -F ',' '{ total += $var } END { printf "%.3e", total/(NR-1) }' $DIR/GOV2_1M-10M_comparison.csv >>$DIR/GOV2_averages.csv
  done
  printf "\n" >>$DIR/GOV2_averages.csv

  printf "'GOV2_AVG_100K-1M'" >>$DIR/GOV2_averages.csv
  for i in $(seq 2 "$NUM_COLUMNS"); do
    printf "," >>$DIR/GOV2_averages.csv
    awk -v var="$i" -F ',' '{ total += $var } END { printf "%.3e", total/(NR-1) }' $DIR/GOV2_100K-1M_comparison.csv >>$DIR/GOV2_averages.csv
  done
  printf "\n" >>$DIR/GOV2_averages.csv

  cat $DIR/DNA_5GRAM_URL_comparison.csv >$DIR/comparison.csv
  cat $DIR/GOV2_averages.csv >>$DIR/comparison.csv
  sed -n '6p' <$DIR/GOV2_10M-_comparison.csv >>$DIR/comparison.csv
}
