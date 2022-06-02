#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=pheno_binarize
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=4G
#SBATCH --chdir /data5/deepro/ukbiobank/preprocess/rarecomb_pheno_prepare/src
#SBATCH -o /data5/deepro/ukbiobank/preprocess/rarecomb_pheno_prepare/slurm/logs/out_download_%a.log
#SBATCH -e /data5/deepro/ukbiobank/preprocess/rarecomb_pheno_prepare/slurm/logs/err_download_%a.log
#SBATCH --array 1-4

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data5/deepro/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data5/deepro/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/data5/deepro/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/data5/deepro/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate ukbiobank

echo `date` starting job on $HOSTNAME
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data5/deepro/ukbiobank/preprocess/rarecomb_pheno_prepare/slurm/files/binarize_pheno_types.txt)

echo $LINE
python /data5/deepro/ukbiobank/preprocess/rarecomb_pheno_prepare/src/0_binarize_phenos.py $LINE 

echo `date` ending job
