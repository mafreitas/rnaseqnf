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
    path 'genome.fa' from params.genome_fasta 
    path 'genes.gtf' from params.genome_gtf

    output:
    path 'star_iundex' into star_index_ch

    script:
    """
    mkdir star-index 
    STAR --runMode genomeGenerate --runThreadN $task.cpus \
        --sjdbOverhang 99 \
        --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf --genomeDir star_iundex/

    """
}
star_index_ch_value = star_index_ch.first()

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
    publishDir "${params.outdir}/trim_galore/${sample_id}", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(reads) from read_pairs3_ch

    output:
    set val("$sample_id"), file("*_trimmed*.fq.gz") into trimgalore_reads_salmon, trimgalore_reads_star
    file "*trimming_report.txt" into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

    script:
    """
    trim_galore --paired --fastqc --gzip --basename ${sample_id}_trimmed ${reads}
    ls
    """
}

process salmon_quant {
    publishDir "$params.outdir/salmon/", mode:'copy'

    cpus 6

    input:
    path index from index_ch_value
    tuple val(sample_id), path(reads) from trimgalore_reads_salmon

    output:
    path(sample_id) into quant_ch

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}

process star_align {
    publishDir "$params.outdir/star_aligned/${sample_id}", mode:'copy'

    cpus 6

    input:
    path index from star_index_ch
    path 'genes.gtf' from params.genome_gtf
    tuple val(sample_id), path(reads) from trimgalore_reads_star

    output:
    set val("$sample_id"), file("*Aligned.out.bam") into star_aligned_bam

    script:
    """
    STAR \
        --genomeDir ${index} \
        --sjdbGTFfile genes.gtf \
        --runThreadN ${task.cpus} \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outFileNamePrefix ${sample_id}_ \
        --readFilesIn ${reads[0]} ${reads[1]} \
        --readFilesCommand zcat \
        --outSAMtype BAM Unsorted \
        --outReadsUnmapped Fastx \
        --outSAMmode Full
    """
}


process feature_count {
    publishDir "$params.outdir/feature_count/${sample_id}", mode:'copy'

    cpus 6

    input:
    path 'genes.gtf' from params.genome_gtf
    tuple val(sample_id), path(bam) from star_aligned_bam

    output:
    set val("$sample_id"), file("${sample_id}_counts.csv") into feature_counts

    script:
    """
    featureCounts \
        -T ${task.cpus} \
        -t exon \
        -g gene_id \
        -a genes.gtf \
        -o ${sample_id}_counts.csv \
        ${bam}
    """
}

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