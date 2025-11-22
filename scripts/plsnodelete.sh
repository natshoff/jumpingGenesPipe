#!/bin/bash
#SBATCH -n 4
#SBATCH -p main
#SBATCH -N 1
#SBATCH --mem=50G

# NP Hofford Mar 3 2024
# This script runs the plsnodelete bash_rc call on the directory you choose so things are not deleted
# Use for scratch and call using sbatch plsnodelete.sh /pathtodirectory

cd $1

find . -exec touch {} \;

echo 'This is so great how now everything isnt going to get deleted. I love that.'
