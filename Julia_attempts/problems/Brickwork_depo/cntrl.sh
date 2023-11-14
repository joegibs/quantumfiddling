#!/bin/bash

for n in {1..5}
do
    sbatch sbfile.sh
    sleep 0.1
done