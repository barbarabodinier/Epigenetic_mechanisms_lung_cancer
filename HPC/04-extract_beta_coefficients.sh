#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=30gb
#PBS -N extracting_beta_80
#PBS -q med-bio
#PBS -J 1-2

cd /rds/general/user/bbodinie/ephemeral/methylation_lung_cancer/Scripts
module load anaconda3/personal

ttd=0
Rscript extract_beta_coefficients.R $PBS_ARRAY_INDEX $ttd
