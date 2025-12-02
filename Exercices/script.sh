#!/bin/bash

FILE="ex5_red"
for thread in $(seq 1 1 16)
do
    export OMP_NUM_THREADS=$thread
    time ./$FILE
done