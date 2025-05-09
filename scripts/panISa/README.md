# ISScreener
Runs insertion sequence detection workflow on Illumina WGS data

# Install
git clone https://github.com/kylegontjes/ISScreener.git

# Usage
module load singularity

module load snakemake

# Dry run
snakemake -s ISScreener.smk --dryrun -p

# Sample command
## Raw fastq files (i.e., trim first)
snakemake -s ISScreener.smk --use-conda --use-singularity -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes}  -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem}  --output=slurm_out/slurm-%j.out" --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000

## Pre-trimmed fastq files  
snakemake -s ISScreener_pretrimmed.smk --use-conda --use-singularity -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes}  -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem}  --output=slurm_out/slurm-%j.out" --conda-frontend conda --cluster-config config/cluster.json --configfile config/config_pretrimmed.yaml --latency-wait 1000

# To run many isolates at the same time, consider using the sbat file
sbatch ISScreener.sbat 

# Curating the samples_list.csv file
## Given a path to a directory with files format

path="/nfs/turbo/umms-esnitkin/Project_Penn_KPC/Sequence_data/illumina_fastq/clean_fastq_qc_pass_samples"

sample_id="sample_id"

sample_names=$(ls -1 $path | grep _R1 |  cut -d. -f1 | sed 's/_R1//' | sed 's/_trim_paired//')

echo -e\n $sample_id $sample_names | tr ' ' '\n' > config/sample.tsv