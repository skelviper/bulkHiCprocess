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

"""
SOFTWARES AND DATA
note that other softwares like cutadapt BWA etc. is also needed, recommend install them with conda
"""

#mm10
BWA_REF_GENOME="/share/home/zliu/share/Data/public/ref_genome/mouse_ref/M23/bwamem2/genome.fa",

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
    output:
        bamfile = "processed/mapping/{sample}.aln.bam"
    threads: 20
    resources:
        nodes = 20
    params:
        extra=r"-R '@RG\tID:{sample}\tPL:ILLUMINA\tSM:{sample}'",
    shell:"""
        set +u
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate py3
        set -u

        bwa-mem2 mem -5SP -t{threads} {params.extra} {BWA_REF_GENOME} {input.DNA1} {input.DNA2} | samtools sort -@{threads} -n -o {output.bamfile} -

        set +u
        conda deactivate
        set -u
        """

rule bam2pairs:
    input:
    	bamfile = rules.bwa_map.output.bamfile,
    	chromsizes = "/share/home/zliu/share/Data/public/ref_genome/mouse_ref/M23/raw_data/chr.len",
    output:
    	dedupPairs = "processed/pairs/{sample}.pairs.gz",
    threads: 15
    resources:
        nodes = 15
    shell:"""
    	set +u
        source activate
        conda activate hic
        set -u

        pairtools parse --chroms-path {input.chromsizes} {input.bamfile} | \
        pairtools sort --nproc 10 --tmpdir=./ | pairtools dedup | pairtools split --output-pairs {output.dedupPairs}

        set +u
        conda deactivate
        set -u
    """
