#!/bin/bash
#SBATCH -J methylation_R
#SBATCH -p rna          
#SBATCH -N 1                 
#SBATCH --mem=200G                 
#SBATCH --cpus-per-task=30   
#SBATCH -t 5-00:00:00        
#SBATCH -o methylation_R.out     
#SBATCH --mail-type=ALL
#SBATCH --mail-user=24211510018@m.fudan.edu.cn


Rpath="/home/u24211510018/.conda/envs/Rlibrary/bin/Rscript"



$Rpath ./Decon.R


# python_path="/home/u24211510018/.conda/envs/DL/bin/python"


# $python_path ./Decon.py