#!/bin/bash

cd 1D/BiKappa/
sbatch slurm_BiKappa_TestScript.slurm
cd ../BiMaxwellian/
sbatch slurm_BiMax_p6_TestScript.slurm
cd ../Kappa/
sbatch slurm_Kappa_p3_TestScript.slurm
cd ../Maxwellian/
sbatch slurm_Max_TestScript.slurm
cd ../IncompleteMax/
sbatch slurm_IMax_TestScript.slurm
