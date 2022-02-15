#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=70gb
#PBS -N lasso
#PBS -q med-bio
#PBS -J 1-50

cd /rds/general/user/bbodinie/ephemeral/methylation_lung_cancer/Scripts
module load anaconda3/personal

Rscript stability.R {model_id_input} $PBS_ARRAY_INDEX {ttd_input}
