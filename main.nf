//Genome Information
params.genomes = "/fs/project/PAS1353/genomes"
params.transcriptome = "${params.genomes}/gencode.v35.transcripts.fa.gz"

params.genome_version = "${params.genomes}/Homo_sapiens/NCBI/GRCh38"
params.genome_fasta = "${params.genome_version}/Sequence/WholeGenomeFasta/genome.fa"
params.genome_gtf = "${params.genome_version}/Annotation/Genes/genes.gtf"

params.reads = "/fs/project/PAS1353/rna_seq_Tsichlis/FASTQ_files/Sample*_R{1,2}_*{1,2}.fq.gz"

params.outdir = "results"
params.multiqc = "multiqc"

Channel
    .fromFilePairs( params.reads, checkExists:true )
    .into { read_pairs_ch; read_pairs2_ch; read_pairs3_ch }

Channel
    .fromPath( params.transcriptome )
    .set { transcriptome }



process salmon_index {
    cpus 26
    input:
    path 'transcriptome.fa' from transcriptome

    output:
    path 'index' into index_ch

    script:
    """
    salmon index --keepDuplicates --threads $task.cpus -t transcriptome.fa -i index --gencode
    """
}
index_ch_value = index_ch.first()

process star_index {
    cpus 26
    input:
    path 'genome.fa' from genome_fasta
    path 'genes.gtf' from genome_gtf

    output:
    path 'star_iundex' into star_iundex_ch

    script:
    """
    mkdir star-index 
    STAR --runMode genomeGenerate --runThreadN $task.cpus \
        --sjdbOverhang 99 \
        --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf --genomeDir star_iundex/

    """
star_iundex_ch_value = star_iundex_ch_value.first()

process fastqc {
     publishDir "$params.outdir/fastqc/$sample_id", mode:'copy'
     cpus 4
     
     input:
     tuple val(sample_id), path(reads) from read_pairs2_ch

     output:
     path "fastqc_${sample_id}_logs" into fastqc_ch

     script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process trim_galore {
    publishDir "${params.outdir}/trim_galore/${name}", mode: 'copy'
    cpus 4

    input:
    tuple val(name), path(reads) from read_pairs3_ch

    output:
    sset val(name), file("*fq.gz") into trimgalore_reads_salmon, trimgalore_reads_star
    file "*trimming_report.txt" into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

    script:
    """
    trim_galore --paired --fastqc --gzip ${reads}
    ls
    """
}

process quant {
    publishDir "$params.outdir/salmon/", mode:'copy'

    cpus 6

    input:
    path index from index_ch_value
    tuple val(pair_id), path(reads) from trimgalore_reads_salmon

    output:
    path(pair_id) into quant_ch

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}



// NOTE TO CHANGE INPUT FILE

// process multiqc {
//     publishDir "$params.outdir/multiqc/", mode:'copy'
    
//     input:
//     path 'data*/*' from quant_ch.mix(fastqc_ch).collect()

//     output:
//     path 'multiqc_report.html'

//     script:
//     """
//     multiqc -v .
//     """
// }




              
// rule trimming:
//     input:
//         fq = "fastqs/{datafile}.fastq.gz"
//     output:
//         fq = "trimmed_fastqs/{datafile}_trimmed.fq.gz"
//     singularity:
//         "docker://mfreitas/rnaseq_star:latest"
//     threads:
//         1
//     shell:
//         "trim_galore "
//         "-q 20 "
//         "-o trimmed_fastqs "
//         "{input.fq}"

// rule ptrimming:
//     input:
//         fq1 = "../parthun_001_rnaseq/paired_fastqs/{datafile}_R1_001.fastq.gz",
//         fq2 = "../parthun_001_rnaseq/paired_fastqs/{datafile}_R2_001.fastq.gz",       
//     output:
//         fq1 = "trimmed_paired_fastqs/{datafile}_R1_001_val_1.fq.gz",
//         fq2 = "trimmed_paired_fastqs/{datafile}_R2_001_val_2.fq.gz"
//     singularity:
//         "docker://mfreitas/rnaseq_star:latest"
//     shell:
//         "trim_galore "
//         "-q 20 "
//         "--paired "
//         "-o trimmed_paired_fastqs "
//         "{input.fq1} {input.fq2}"

// rule fastqc:
//     input:
//         fq = "fastqs/{datafile}.fastq.gz"
//     output:
//         fq = "fastqcout/{datafile}_fastqc.html"
//     singularity:
//         "docker://mfreitas/rnaseq_star:latest"
//     shell:
//         "fastqc "
//         "{input.fq} "
//         "-o fastqcout"
        
// rule pfastqc:
//     input:
//         fq = "../parthun_001_rnaseq/paired_fastqs/{datafile}.fastq.gz",
//     output:
//         fq = "paired_fastqcout/{datafile}_fastqc.html"
//     singularity:
//         "docker://mfreitas/rnaseq_star:latest"
//     shell:
//         "fastqc "
//         "{input.fq} "
//         "-o paired_fastqcout"
      
// rule tfastqc:
//     input:
//         fq = "trimmed_fastqs/{datafile}_trimmed.fq.gz"
//     output:
//         fq = "trimmed_fastqcout/{datafile}_trimmed_fastqc.html"
//     singularity:
//         "docker://mfreitas/rnaseq_docker://latest"
//     shell:
//         "fastqc "
//         "{input.fq} "
//         "-o trimmed_fastqcout"
     
// rule tpfastqc:
//     input:
//         fq1 = "trimmed_paired_fastqs/{datafile}_R1_001_val_1.fq.gz",
//         fq2 = "trimmed_paired_fastqs/{datafile}_R2_001_val_2.fq.gz"
//     output:
//         fq1 = "trimmed_paired_fastqcout/{datafile}_R1_001_val_1_fastqc.html",
//         fq2 = "trimmed_paired_fastqcout/{datafile}_R2_001_val_2_fastqc.html"
//     singularity:
//         "docker://mfreitas/rnaseq_star:latest"
//     shell:
//         "fastqc "
//         "{input.fq1} {input.fq2} "
//         "-o trimmed_paired_fastqcout"  
           
// rule align:
//     input:
//         fq = "trimmed_fastqs/{datafile}_trimmed.fq.gz",
//         genome = "reference/Genome"
//     output:
//         fq = "star_output/{datafile}_Aligned.out.bam"
//     threads:
//         8
//     singularity:
//         "docker://mfreitas/rnaseq_star:latest"
//     shell:
//         "STAR "
//         "--genomeDir reference "
//         "--sjdbGTFfile {genome_gtf} "
//         "--runThreadN {threads} "
//         "--outSAMstrandField intronMotif "
//         "--outFilterIntronMotifs RemoveNoncanonical "
//         "--outFileNamePrefix star_output/{wildcards.datafile}_ "
//         "--readFilesIn {input.fq} "
//         "--readFilesCommand zcat "
//         "--outSAMtype BAM Unsorted "
//         "--outReadsUnmapped Fastx "
//         "--outSAMmode Full"

// rule palign:
//     input:
//         fq1 = "trimmed_paired_fastqs/{datafile}_R1_001_val_1.fq.gz",
//         fq2 = "trimmed_paired_fastqs/{datafile}_R2_001_val_2.fq.gz",
//         genome = "reference/Genome"
//     output:
//         fq = "star_output/{datafile}_Aligned.out.bam"
//     threads:
//         8
//     singularity:
//         "docker://mfreitas/rnaseq_star:latest"
//     shell:
//         "STAR "
//         "--genomeDir reference "
//         "--sjdbGTFfile {genome_gtf} "
//         "--runThreadN {threads} "
//         "--outSAMstrandField intronMotif "
//         "--outFilterIntronMotifs RemoveNoncanonical "
//         "--outFileNamePrefix star_output/{wildcards.datafile}_ "
//         "--readFilesIn {input.fq1} {input.fq2} "
//         "--readFilesCommand zcat "
//         "--outSAMtype BAM Unsorted "
//         "--outReadsUnmapped Fastx "
//         "--outSAMmode Full"

// rule featurecount:
//     input:        
//         fq = "star_output/{datafile}_Aligned.out.bam"
//     output:
//         csv = "counts/{datafile}_counts.csv"
//     threads:
//         8
//     singularity:
//         "docker://mfreitas/rnaseq_star:latest"
//     shell:
//         "featureCounts "
//         "-T {threads} "
//         "-t exon "
//         "-g gene_id "
//         "-a {genome_gtf} "
//         "-o {output.csv} "
//         "{input.fq}" 
        
// rule featurecountall:
//     input:        
//         fqs = expand("star_output/{sample}_Aligned.out.bam",sample=SAMPLES+PESAMPLES)
//     output:
//         csv = "counts/all_counts.csv"
//     threads:
//         8
//     singularity:
//         "docker://mfreitas/rnaseq_star:latest"
//     shell:
//         "featureCounts "
//         "-T {threads} "
//         "-t exon "
//         "-g gene_id "
//         "-a {genome_gtf} "
//         "-o {output.csv} "
//         "{input.fqs}" 
        
// rule multiqc:
//     input:        
//         "counts/all_counts.csv"
//     output:
//         "multiqc_report.html"
//     threads:
//         8
//     singularity:
//         "docker://mfreitas/rnaseq_star:latest"
//     shell:
//         "multiqc ." 
