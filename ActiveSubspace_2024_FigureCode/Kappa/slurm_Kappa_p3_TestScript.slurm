#!/bin/bash -x

# NOTE: all these #SBATCH lines are not commented (though this line is, as is anything with ##SBATCH)

#SBATCH --job-name=“Kappa_Sampling_Error_Test”
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

##SBATCH —-account=181218142210
#SBATCH --exclusive # this is for an exclusive run on one of the 'pankavic' nodes
#SBATCH --export=ALL
#SBATCH --time=5:00:00 
# this is the maximum time limit for a job (DD-HH:MM:SS)
                           # any longer than 6 days and the scheduler will kill
                           # it at the end of that time period
			   ##SBATCH --time=05-23:00:00 is typical limit

##SBATCH --mail-type=ALL
##SBATCH --mail-user=email@mines.edu # if you uncomment this, you will get emails when the job begins and ends
#SBATCH -o lastrunp6.out # these two are the out and error files for the batch scheduler
#SBATCH -e lastrunp6.err
##SBATCH -v # this will give more verbose output from the batch scheduler
cd $SLURM_SUBMIT_DIR

# this loads the matlab module--can also just put this in your .bashrc to have it done when you log in

module load apps/matlab/2021b

# run the code and direct all matlab (i.e., command window) output to outfilep6.txt

matlab -nodisplay -nosplash -r TestScript_Kap
