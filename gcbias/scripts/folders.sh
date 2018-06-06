#!/bin/bash

base_path="/paula/2018/kostka"

data_path=$base_path"/data"
mkdir $data_path
cd $data_path
mkdir ./config
mkdir ./input
mkdir ./preparation
cd ./preparation
mkdir ./elements_bed
mkdir ./elements_msa
mkdir ./elements_msa_ss
mkdir ./flanking_region_bed
mkdir ./flanking_region_msa
mkdir ./flanking_region_msa_ss

mkdir ../results

cd ./elements_bed
mkdir ./1_chr
mkdir ./2_elements
cd ../flanking_region_bed
mkdir ./1_chr
mkdir ./2_elements
mkdir ./3_subtract_exclusion
mkdir ./4_subtract_element


mkdir $base_path"/scripts"

