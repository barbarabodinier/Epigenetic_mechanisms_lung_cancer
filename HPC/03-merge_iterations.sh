#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=1:mem=60gb
#PBS -N merging
#PBS -q med-bio
#PBS -J 1-2

cd /rds/general/user/bbodinie/ephemeral/methylation_lung_cancer/Scripts
module load anaconda3/personal

ttd=0
echo $ttd
Rscript merge_iterations.R $PBS_ARRAY_INDEX $ttd
