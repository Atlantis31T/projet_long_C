#!/bin/bash

mpirun -np $1 clusters0 $2 | tee output
ls cluster.final.* | wc -l > nbclusters

