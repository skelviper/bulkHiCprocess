####################################
#    Fastq2pairs for bulk Hi-C     #
#@author Z Liu                     #
#@Ver 0.1                          #
#@date 2020/6/09                   #
####################################

#############CONFIG#################

import os

def listdir_nohidden(path):
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f.split(sep='_')[0]
            
SAMPLES = list(set(listdir_nohidden("./Rawdata")))
#SAMPLES = [""]

configfile: "bulkHiCprocess/config.yaml"

#############END_CONFIG#############

#############RULE_ALL###############
"""
decide what you need for your down stream analysis.
"""
rule all:
    input:
    	expand("processed/pairs/{sample}.pairs.gz",sample=SAMPLES),


rule bwa_map:
    input:
        DNA1="Rawdata/{sample}_1.fastq.gz",
        DNA2="Rawdata/{sample}_2.fastq.gz",
        BWA_REF_GENOME=config["refs"][config["ref_genome"]]["BWA_REF_GENOME"],
    output:
        bamfile = "processed/mapping/{sample}.aln.bam"
    threads: 30
    resources:
        nodes = 30
    params:
        extra=r"-R '@RG\tID:{sample}\tPL:ILLUMINA\tSM:{sample}'",
    shell:"""
        set +u
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate py3
        set -u

        bwa-mem2 mem -5SP -t{threads} {params.extra} {input.BWA_REF_GENOME} {input.DNA1} {input.DNA2} | samtools sort -@{threads} -n -o {output.bamfile} -

        set +u
        conda deactivate
        set -u
        """

rule bam2pairs:
    input:
    	bamfile = rules.bwa_map.output.bamfile,
    	chromsizes = config["refs"][config["ref_genome"]]["chr_len"],
    output:
    	dedupPairs = "processed/pairs/{sample}.pairs.gz",
    threads: 30
    resources:
        nodes = 30
    shell:"""
    	set +u
        source activate
        conda activate hic
        set -u

        pairtools parse --chroms-path {input.chromsizes} {input.bamfile} | \
        pairtools sort --nproc 30 --tmpdir=./ | pairtools dedup | pairtools split --output-pairs {output.dedupPairs}

        set +u
        conda deactivate
        set -u
    """
