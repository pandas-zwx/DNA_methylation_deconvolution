#!/bin/bash
#SBATCH -J grid_search
#SBATCH -p cell     
#SBATCH -N 1                 
#SBATCH --mem=80G                 
#SBATCH --cpus-per-task=20   
#SBATCH -t 5-00:00:00        
#SBATCH -o grid_search.out     
#SBATCH --mail-type=ALL
#SBATCH --mail-user=24211510018@m.fudan.edu.cn


~/.conda/envs/translation/bin/python grid_search_for_decon.py