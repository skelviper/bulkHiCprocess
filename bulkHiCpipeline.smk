import os
def listdir_nohidden(path):
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f.replace(".pairs.gz","")

SAMPLES = listdir_nohidden("./pairs")

##########reference file config################
comaprtment_ref_track = "/share/home/zliu/share/Data/public/ref_genome/GRCh37/raw_data/hg19.gccontent.withchr.100k.txt"
###############################################
# call compartment and TAD like Hi-C 3.0
# draw saddle plot and calc comapartment strength



rule all:
    input:
        expand("TADs/{sample}.tad.bed",sample=SAMPLES),
        expand("compartment/{sample}.compartment.100k.cis.vecs.tsv",sample=SAMPLES),
        expand("expected/{sample}.100k.expected.tsv",sample=SAMPLES),
        expand("TADres/{sample}.TADres.tsv",sample=SAMPLES),
        expand("saddle/{sample}.saddle.saddledump.npz",sample=SAMPLES),
        "saddle/compartmentStrength.txt",  

rule pairs2cool:
    input:
        pairs = "pairs/{sample}.pairs.gz",
    output:
        balancedCool = "cools/{sample}.5k.balanced.cool",
    params:
        resolution = 5000,
    threads: 8
    shell:"""
        set +u
        source activate
        conda activate hic
        set -u

        cooler cload pairs -c1 2 -c2 4 -p1 3 -p2 5 /share/Data/public/ref_genome/GRCh37/raw_data/hg19.chr.withchr.len:{params.resolution} {input.pairs} {output.balancedCool}
        cooler balance {output.balancedCool}

        set +u
        conda deactivate
        set -u
    """

rule generate_mcool:
    input:
        balancedCool = rules.pairs2cool.output.balancedCool,
    output:
        mcool = "mcools/{sample}.balanced.mcool"
    threads: 10
    shell:"""
    
        set +u
        source activate
        conda activate hic
        set -u
        
        cooler zoomify -p {threads} {input.balancedCool} -r 5000,10000,20000,40000,100000,200000,500000,1000000 --balance -o {output.mcool} 

        set +u
        conda deactivate
        set -u
    """
    
rule call_compartment:
    input:
        rules.generate_mcool.output.mcool,
    output:
        compartment = "compartment/{sample}.compartment.100k.cis.vecs.tsv"
    params:
        resolution=100000
    shell:"""
        set +u
        source activate
        conda activate hic
        set -u

        cooltools call-compartments {input}::/resolutions/{params.resolution} --reference-track {comaprtment_ref_track} -o ./compartment/{wildcards.sample}.compartment.100k --bigwig

        set +u
        conda deactivate
        set -u
    """

rule compute_expected:
    input:
        mcool = rules.generate_mcool.output.mcool
    output:
        #expected5k = "expected/{sample}.5k.expected.tsv"
        expected100k = "expected/{sample}.100k.expected.tsv"

    threads: 10
    params:
        resolution=100000
    shell:"""
        set +u
        source activate
        conda activate hic
        set -u

        cooltools compute-expected {input.mcool}::/resolutions/{params.resolution}  -p {threads} -o {output.expected100k}

        set +u
        conda deactivate
        set -u
    """

rule call_tad:
    input:
        coolpath = rules.pairs2cool.output.balancedCool,
    output:
        insulationScore = "insulation/{sample}.standardTAD.tsv",
    params:
        windowSize = 10000,
    shell:"""
        set +u
        source activate
        conda activate hic
        set -u
        
        cooltools diamond-insulation {input.coolpath} --ignore-diags 2 --window-pixels 40 > {output.insulationScore}
        
        set +u
        conda deactivate
        set -u
    """

rule compartmentVecs2bed:
    input:
        rawCom = rules.call_compartment.output.compartment,
    output:
        compartmentbed = "compartment/{sample}.compartment.100k.bed",
        tmpfile = temp("compartment/{sample}.temp.bed")
    shell:"""
        set +u
        source activate
        conda activate R
        set -u

        Rscript compartmentVecs2bed.R {input.rawCom} {output.tmpfile}

        bedtools merge -s -d 3 -i {output.tmpfile} -o distinct -c 5 > {output.compartmentbed}

        set +u
        conda deactivate
        set -u
    """

rule cleanTADs:
    input:
        rawTADs = rules.call_tad.output.insulationScore,
        comp = rules.compartmentVecs2bed.output.compartmentbed,
    output:
        cleanedTADs = "TADs/{sample}.tad.bed"

    shell:"""
        set +u
        source activate
        conda activate R
        set -u

        Rscript cleanTADs.R {input.comp} {input.rawTADs} {output.cleanedTADs}

        set +u
        conda deactivate
        set -u
    """

rule compute_saddle:
    input:
        expected = rules.compute_expected.output.expected100k,
        mcool = rules.generate_mcool.output.mcool,
        compartment = "compartment/{sample}.compartment.100k.cis.vecs.tsv",
    output:
        saddleplot = "saddle/{sample}.saddle.png",
        saddle_strength = "saddle/{sample}.saddle.saddledump.npz",
    params:
        resolution = 100000
    shell:"""
        set +u
        source activate
        conda activate hic
        set -u
        cooltools compute-saddle {input.mcool}::/resolutions/{params.resolution} \
        ./compartment/{wildcards.sample}.compartment.100k.cis.vecs.tsv::E1 {input.expected} \
        -o ./saddle/{wildcards.sample}.saddle --fig png --strength --qrange 0.025 0.975
        set +u
        conda deactivate
        set -u
    """

rule compute_strength:
    input:
        expand(rules.compute_saddle.output.saddle_strength,sample = SAMPLES)
    output:
        "saddle/compartmentStrength.txt"
    run:
        with open(output[0],'w') as out:
            for file in input:
                out.write(file)
                out.write(np.load("".join([celltype,".saddle.saddledump.npz"]))['saddle_strength'][10])
        print("All done!")

