#!/bin/bash
#

#SBATCH --job-name=sampling
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jacqueline.rudolph@jhu.edu
#SBATCH --output=./results/poisson_option.out

#SBATCH --partition=lau
#SBATCH --nodelist=compute-135
#SBATCH --time=7-0:00:00
#SBATCH --mem=memory
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --verbose

t=$(date)
echo -e "Start time: $t \n"
d=$(pwd)
echo -e "Working directory: $d \n"

echo -e "Analysis: option \n"

module load R
#Rscript ./code/samp_2_poisson-divide.R option
Rscript ./code/samp_2_poisson-divide-repeat.R option

echo
sstat -a -o JobID,MaxVMSizeNode,MaxVMSize,AveVMSize,MaxRSS,AveRSS,MaxDiskRead,MaxDiskWrite,AveCPUFreq,TRESUsageInMax -j ${SLURM_JOB_ID}
echo
t=$(date)
echo "End time: $t"
