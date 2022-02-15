#PBS -l walltime=3:00:00
#PBS -l select=1:ncpus=1:mem=30gb
#PBS -N lasso_cv
#PBS -q med-bio
#PBS -J 1-2

cd /rds/general/user/bbodinie/ephemeral/methylation_lung_cancer/Scripts
module load anaconda3/personal

Rscript cv_lasso.R $PBS_ARRAY_INDEX
