#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "data/*.{1,2}.fastq.gz" // it is assumed that paired-end sequencing technology is being used.
params.genome = "references/human_genome.fa"
params.vep_dir = "/home/user926/.vep"
params.threads = 2

def checkExecutable(exec_name, hint, alt_path = null, version_cmd = null) {
    def proc = "which ${exec_name}".execute()
    proc.waitFor()
    if (proc.exitValue() != 0) {
        println "[FAIL] '${exec_name}' not found in PATH"
        println "       ➤ Hint: ${hint}"
        if (alt_path != null) {
            println "       ➤ If already installed locally, add to PATH using:"
            println "         export PATH=\"${alt_path}:\$PATH\""
        }
        println ""
        return false
    } else {
        String version_info = ""
        if (version_cmd) {
            try {
                def version_proc = ["bash", "-c", version_cmd].execute()
                version_proc.waitFor()
                version_info = version_proc.in.text.readLines().find { it } ?: "(no output)"
                version_info = version_info.trim()
            } catch (Exception e) {
                version_info = "(version unavailable)"
            }
        }
        println "[ OK ] '${exec_name}' is available - ${version_info}"
        return true
    }
}

workflow {
    main:

    def all_ok = true

    def tools = [
        ['fastqc',     "conda install -c bioconda fastqc",     null, "fastqc --version"],
        ['fastp',      "conda install -c bioconda fastp",      null, "fastp --version"],
        ['bwa',        "conda install -c bioconda bwa",        null, "bwa 2>&1 | grep -i version | head -n 1"],
        ['samtools',   "conda install -c bioconda samtools",   null, "samtools --version | head -n 1"],
        ['gatk',       "conda install -c bioconda gatk4",      null, "gatk --version"],
        ['vcf2maf.pl', "https://github.com/mskcc/vcf2maf",     "$HOME/vcf2maf", "vcf2maf.pl --help | head -n 1"]
    ]

    for (tool in tools) {
        def (name, install_cmd, alt_path, version_cmd) = tool
        if (!checkExecutable(name, install_cmd, alt_path, version_cmd)) {
            all_ok = false
        }
    }

    if (!file(params.genome).exists()) {
        println "[FAIL] Genome reference not found: ${params.genome}"
        println "       ➤ Hint: Place your genome FASTA file at the specified path."
        println "         curl -O ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        println "         gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        println "         mv Homo_sapiens.GRCh38.dna.primary_assembly.fa ${params.genome}"
        println "         samtools faidx ${params.genome}"
        println "         gatk CreateSequenceDictionary -R ${params.genome}"
        all_ok = false
    } else {
        println "[ OK ] Genome reference found: ${params.genome}"
    }

    if (!file(params.vep_dir).exists()) {
        println "[FAIL] VEP directory not found: ${params.vep_dir}"
        println "       ➤ Hint: Install and configure VEP, then set params.vep_dir to its path."
        all_ok = false
    } else {
        println "[ OK ] VEP directory found: ${params.vep_dir}"
    }

    Channel.fromFilePairs(params.reads, flat: true)
           .collect()
           .view { pairs ->
                if (pairs.isEmpty()) {
                    println "[FAIL] No reads found matching: ${params.reads}"
                    println "       ➤ Hint: Make sure your FASTQ files are in the 'data/' directory and follow the *_1.fastq.gz, *_2.fastq.gz pattern."
                    all_ok = false
                } else {
                    println "[ OK ] Found ${pairs.size()} read pair(s) matching: ${params.reads}"
                }
            }

    if (!all_ok) {
        println "\nOne or more checks failed. Please review the hints above."
        exit 1
    } else {
        println "\nAll checks passed successfully!"
    }
}
