params.reads = "/fs/project/PAS1353/rna_seq_Tsichlis/FASTQ_files/Sample*_R{1,2}_*{1,2}.fq.gz"
params.transcriptome = "/fs/project/PAS1353/genomes/gencode.v35.transcripts.fa.gz"
params.outdir = "results"
params.multiqc = "multiqc"


log.info """\
 R N A S E Q - N F   P I P E L I N E
 ===================================
 transcriptome: ${params.transcriptome}
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """


Channel
    .fromFilePairs( params.reads, checkExists:true )
    .into { read_pairs_ch; read_pairs2_ch; read_pairs3_ch }

Channel
    .fromPath( params.transcriptome )
    .set { transcriptome }


process index {
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


//read_pairs_ch.view()
index_ch_value = index_ch.first()


//read_pairs_ch.view()
process quant {
    publishDir "$params.outdir/salmon/", mode:'copy'

    cpus 6

    input:
    path index from index_ch_value
    tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    path(pair_id) into quant_ch

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}


 process fastqc {
     publishDir "$params.outdir/fastqc_1/$sample_id", mode:'copy'
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

process trim_galore_1 {
    publishDir "${params.outdir}/trim_galore_1/trim_1", mode: 'copy'
    cpus 4

    input:
    tuple val(pair_id), path(reads) from read_pairs3_ch

    output:
    set val(name), file("*fq.gz") into  trimgalore_reads
    file "*trimming_report.txt" into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

    script:
    """
    trim_galore --paired --fastqc --gzip ${reads}
    ls
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
