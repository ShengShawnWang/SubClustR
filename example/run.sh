#!/bin/bash

Rscript scripts/SubClustR.R \
  --expression_data data/example_expression.csv \
  --genes "GENE1,GENE2,GENE3" \
  --patience 25 \
  --stabilization_threshold 0.01 \
  --fdr 0.01 \
  --out example/



Rscript /projects/b1080/sw/project/personal_github/scripts/SubClustR.R \
  --expression_data /projects/b1080/sw/project/personal_github/example/raw_pc_data.csv \
  --genes all \
  --patience 20 \
  --stabilization_threshold 0.005 \
  --fdr 0.01 \
  --out /projects/b1080/sw/project/personal_github/example/


