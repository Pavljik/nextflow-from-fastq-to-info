#!/usr/bin/env nextflow
nextflow.enable.dsl=2 // It enables DSL2, the modern version of the Nextflow language
                      // I don't really understand why the developers implemented the new version this way
                      // perhaps it's necessary to support legacy pipelines
                      // Nevertheless, this line is important for the code to work correctly

params.sample_list = "samples.txt"
params.reads     = "data/*.{1,2}.fastq.gz"
params.genome    = "${baseDir}/references/genome.fa"
params.vep_dir   = "/home/user/.vep" // Here you need to specify the path to VEP; otherwise, it doesn't work
params.vep_path  = "/home/user/vep/ensembl-vep-release-111" // also
params.threads   = 4
params.outdir    = "results"
params.ref_name  = "GRCh38"


workflow {
                                      // Here we create
    sample_ids = Channel              // a data stream that will "flow through" the pipeline :D
        .fromPath(params.sample_list) // open a text file
        .splitText()                  // split it into line
        .map { it.trim() }            // remove special characters (such as spaces, tabs, and line breaks)
        .toSet()                      // convert the lines into a `set` of elements

    Channel 
        .fromFilePairs(params.reads, size: 2) // the  method creates a channel emitting the file pairs 
                                              // matching a glob pattern provided by the user
                                              // this is convenient if we work with paired fastq from Illumina
        .filter { sample_id, _ -> sample_ids.contains(sample_id) } // take only those samples that are in the file samples.txt
        .set { read_pairs }

    read_pairs.view { pair -> "--> ${pair[0]} --> ${pair[1]*.name.join(', ')}" } // print()

    qc_fastqc(read_pairs)             // 

    // Process chain
    trimmed = fastp(read_pairs)
    aligned = align_bwa(trimmed)
    sorted  = sort_bam(aligned)
    dedup   = mark_duplicates(sorted) //
    vcf     = call_variants(dedup)
    vep_out = run_vep(vcf)
    maf     = convert_to_maf(vep_out) //
    filtered= filter_maf(maf)
}

process qc_fastqc {
    tag "$sample_id"
    input:
    tuple val(sample_id), file(reads)

    output:
    file("*_fastqc.zip")

    // a Nextflow command that controls where the process output files are placed
    // Nextflow copies files from `work/` (cache)
    publishDir "${params.outdir}/fastqc", mode: 'copy' // `link`, `move`, `copyNoFollow`

    script:
    // the simplest utility for assessing the quality of FASTQ files ;)
    // -t --> how many threads to use
    """
    fastqc -t ${params.threads} ${reads[0]} ${reads[1]}
    """
}

process fastp {
    tag "$sample_id"
    input:
    tuple val(sample_id), file(reads)

    output:
    tuple val(sample_id), file("${sample_id}_R1.clean.fastq.gz"), file("${sample_id}_R2.clean.fastq.gz")

    script:
    // -i --> path to the input file R1
    // -I --> path to the Input file R2
    // -o --> path to the output cleaned file R1
    // -O --> path to the Output cleaned file R2
    // -w --> how many threads work
    """
    fastp -i ${reads[0]} -I ${reads[1]} \\
          -o ${sample_id}_R1.clean.fastq.gz \\
          -O ${sample_id}_R2.clean.fastq.gz \\
          -w ${params.threads}
    """
}

process align_bwa {
    tag "$sample_id"
    input:
    tuple val(sample_id), file(r1), file(r2)

    output:
    tuple val(sample_id), file("${sample_id}.bam")

    script:
    // bwa mem --> alignment algorithm
    // @RG is the start of the read group line; \\t is an escaped tab; ID:${sample_id} is the ID of the reading group;
    // SM:${sample_id} is the sample name; PL:ILLUMINA is the platform (e.g. Illumina)
    // if you use bwa for the first time it may seem confusing but then you will get used to it 
    """
    bwa mem -t ${params.threads} \\
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \\
        ${params.genome} ${r1} ${r2} | \\
        samtools view -Sb - > ${sample_id}.bam
    """
    // `|` pipes the output of bwa mem to samtools view
    // `-S` specifies that the input format is SAM; `b` converts to BAM
}

process sort_bam {
    tag "$sample_id"
    input:
    tuple val(sample_id), file(bam)

    output:
    tuple val(sample_id), file("${sample_id}.sorted.bam")

    // samtools sort --> a command from `samtools` sorts the BAM file by coordinates. This is necessary before indexing the BAM file and calling variants
    // -@ --> also shows how many threads can be used such as `-t` or `-w`
    // -o -->  output name
    script:
    """
    samtools sort -@ ${params.threads} -o ${sample_id}.sorted.bam ${bam}
    """
}

process mark_duplicates {
    tag "$sample_id"
    input:
    tuple val(sample_id), file(bam)
    output:
    tuple val(sample_id), file("${sample_id}.dedup.bam"), file("${sample_id}.dedup.bam.bai")

    // -I --> input
    // -O --> output
    // -M --> duplicate statistics
    // samtools index creates the `.bai` index for the `.bam` file so that it can be quickly read in downstream analysis
    script:
    """
    gatk MarkDuplicates -I ${bam} -O ${sample_id}.dedup.bam -M ${sample_id}.metrics.txt
    samtools index ${sample_id}.dedup.bam
    """
}


process call_variants {
    tag "$sample_id"
    input:
    tuple val(sample_id), file(bam), file(bai)
    output:
    tuple val(sample_id), file("${sample_id}.vcf")

    publishDir "${params.outdir}", mode: 'copy'

    // -R specifies the path to the FASTA reference file. Not Reads!
    // -I --> you know that
    // -O --> you know that
    script:
    """
    gatk HaplotypeCaller -R ${params.genome} -I ${bam} -O ${sample_id}.vcf
    """
}

process run_vep {
    tag "$sample_id"
    input:
    tuple val(sample_id), file(vcf)

    output:
    tuple val(sample_id), file("${sample_id}_vep.vcf")

    publishDir "${params.outdir}", mode: 'copy'

// --cache --> use local annotation cache (it is faster but requires 30-100 Gb)
// --offline --> without access the Internet
// --assembly --> GRCh38, GRCh37, GRCm39 etc.
// --vcf --> output format: VCF with annotation
    script:
    """
    vep --cache --offline \\
        --dir_cache ${params.vep_dir} \\
        --assembly ${params.ref_name} \\
        -i ${vcf} -o ${sample_id}_vep.vcf --vcf
    """
}

process convert_to_maf {
    tag "$sample_id"
    input:
    tuple val(sample_id), file(vep_vcf)

    output:
    file("${sample_id}.maf")

    publishDir "${params.outdir}", mode: 'copy'

    // vcf2maf is a Perl script. To run the script you do not need to know Perl but you need to install it on your laptop
    // BTW Mac OS already has Perl installed. Open a Terminal  and run perl -v to find out which version
    script:
    """
    vcf2maf.pl --input-vcf ${vep_vcf} \\
               --output-maf ${sample_id}.maf \\
               --ref-fasta ${params.genome} \\
               --ncbi-build ${params.ref_name} \\
               --vep-path ${params.vep_path}
    """
}

process filter_maf {
    tag "$sample_id"

    input:
    tuple val(sample_id), file(maf)

    output:
    file("${sample_id}.filtered.maf")

    publishDir "${params.outdir}", mode: 'copy'

    // this is a great feature. We can set dependencies for the filter_maf.py
    conda:
    """
    channels:
      - defaults
      - conda-forge
    dependencies:
      - python=3.10
      - pandas=2.3.1
      - numpy=2.3.1
      - requests=2.32.4
      - tqdm=4.67.1
      - myvariant=1.0.0
    """

    script:
    """
    python3 filter_maf.py ${sample_id}
    """
}
