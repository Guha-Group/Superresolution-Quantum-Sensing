#!/bin/bash
### The account associated with our lab
#SBATCH --account=ashoka
### use standard partition to get priority. Otherwise use windfall
#SBATCH --partition=standard
### Time allowed for job.
#SBATCH --time=00-48:00:00
### Number of nodes assigned to EACH JOB. Not the total number of nodes used.
#SBATCH --nodes=1
### number of tasks
#SBATCH --ntasks=1
### number of cpus per task
#SBATCH --cpus-per-task=94
### the name of the job
#SBATCH --job-name=ColorCenterSensing
### memory per CPU
#SBATCH --mem-per-cpu=5gb
### total memory requested = memory per cpu x cpus per task. Better not to use if you have ntasks and cpus per task set
### SBATCH -mem=470gb
### email the user once job is done
#SBATCH --mail-user=ndeshler@arizona.edu
### setup job array
#SBATCH --array 1-5
### output for a job array
#SBATCH --output=ColorCenter_BrightnessBiasIndex-%a.out

#########################
### Additional flags ####
#########################
### Request GPU resources. Up to four GPUs may be requested on Puma on a single GPU node with --gres=gpu:1, 2, 3, or 4
### SBATCH --gres=gpu:1


####################################
### USEFUL ENVIRONMENT VARIABLES ###
####################################
### https://blog.ronin.cloud/slurm-job-arrays/
### $SLURM_ARRAY_JOB_ID
### $SLURM_ARRAY_TASK_ID
### $SLURM_ARRAY_TASK_MAX
### $SLURM_ARRAY_TASK_MIN


####################
### USEFUL LINKS ###
####################
### https://blog.ronin.cloud/slurm-job-arrays/



module load matlab

matlab matlab -nodisplay -nosplash -batch "BayesianReceiver_MonteCarlo($SLURM_ARRAY_TASK_ID,'94','TypeII_Receivers_v3')" output_{$SLURM_ARRAY_TASK_ID}.txt
