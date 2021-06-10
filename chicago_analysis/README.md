# CHiCAGO_for-Radicl-seq
This repo contains the scripts used to generate the input files and design files for running chicago algorithm on Radicl-seq data.


## Plot_try.R
* This file was used to generate initial plots of detecting an interaction based on distance of the tag location.

## Gen.R
* This file is used to generate the baitmap and rmap file without removing overlapping regions.
* The inputs required are the initial data.bedpe file, chromosome size file and a text file containing genomic location for the RNAs.

## matrix_try.R
* This file generates the chicago input file.
* Inputs required are .bedpe file,chromosome size file,rmap and baitmap files.
* The scripts also removes overlapping regions in the rmap file.

## chicago_try.R
* Script used for running the chicago algorithm.

## Analysis.R
* This contains the initial analysis done on the output from chicago run.

## Please open up these scripts in Rstudio instead of running from the terminal. The input files need to be provided manually. I will update these scripts soon. 
