#! /bin/sh

#SBATCH --job-name "RunPando"
#SBATCH --mem=16000
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --qos hubbioit
#SBATCH --partition hubbioit
#SBATCH --output=log/run_pando.out
#SBATCH --error=log/run_pando.err
#SBATCH --mail-type=END
#SBATCH --mail-user=stephan.fischer@pasteur.fr


# required modules
module load R/4.4.0
RESOLUTION=$1
TISSUE=$2
ALL_STAGES=$3

# pseudobulk 15-20min, 6GB
# metacells 15min-1h, 4GB 
# single-cell 3h15-4h30, 12GB
Rscript 04_run_pando.R $RESOLUTION $TISSUE $ALL_STAGES

exit 0
