#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Mandatory Params
params.sample_table = "nash_samples_SRA.csv"
params.Project = ""
params.mod_FastP = ""
params.Organism = ""
params.Library_Layout = ""
params.Control_name = "Control"
params.Case_name = "Case"

process Dowload_fastqs {
    publishDir "${params.outdir}/Raw_data/", mode:'copy'
    maxForks 1   
    errorStrategy 'retry'
    maxRetries 3
    tag "Download fastqs on $run_acc"

    input:
    val(run_acc)

    output:  
    tuple val(run_acc), path("${run_acc}_{1,2}.fastq.gz"), emit: fastqs 

    script:
    """
    if [ `ls ../../../${params.outdir}Raw_data/${run_acc}*.fastq.gz > /dev/null ; echo \$?` -eq 0 ]; then cp ../../../${params.outdir}Raw_data/${run_acc}*.fastq.gz . ; else fastq-dump --split-files --gzip ${run_acc}; fi
    """   
}
process FASTQC {
    publishDir "${params.outdir}/FASTQC", mode:'copy', pattern: '*.html'      
    tag "FASTQC on $sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:    
    path "*.html"
    path "*zip", emit: multiqc_input

    script:
    """
    fastqc $reads
    """
}
process MULTIQC {    
    publishDir params.outdir, mode:'copy'  

    input:
    file archivos
    val(name)

    output:
    path "*_multiqc_report.html"
    
    script:
    """    
    multiqc $archivos -n '$name'_multiqc_report.html
    """
}
process FASTP {    //
    maxForks 5    
    publishDir "${params.outdir}/Trimmed", mode:'copy', pattern: '*.fastq.gz' 
    tag "FASTP on $sample_id"    
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads
    
    script:
    if (params.Library_Layout == "PAIRED")    
        """
        fastp --thread $task.cpus $params.mod_FastP --in1 ${reads[0]} --in2 ${reads[1]} --length_required 40 --detect_adapter_for_pe --cut_right --cut_right_window_size 4 --cut_right_mean_quality 30 --trim_tail1=5 --trim_tail2=5 --trim_front1=15 --trim_front2=15 --out1 '$sample_id'_trimmed_R1.fastq.gz --out2 '$sample_id'_trimmed_R2.fastq.gz
        """
    else if (params.Library_Layout == "SINGLE")     
        """
        fastp --thread $task.cpus $params.mod_FastP --in1 ${reads[0]} --length_required 40 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 30 --trim_tail1=5 --trim_front1=15 --out1 '$sample_id'_trimmed_R1.fastq.gz
        """
    else
        error "Undefined Library_Layout mode: $params.Library_Layout"
}
process FASTP_short_reads {    //
    maxForks 5    
    publishDir "${params.outdir}/Trimmed", mode:'copy', pattern: '*.fastq.gz' 
    tag "FASTP on $sample_id"    
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads
    
    script:
    if (params.Library_Layout == "PAIRED")    
        """
        fastp --thread $task.cpus $params.mod_FastP --in1 ${reads[0]} --in2 ${reads[1]} --length_required 30 --detect_adapter_for_pe --cut_right --cut_right_window_size 4 --cut_right_mean_quality 30 --trim_front1=10 --trim_front2=10 --out1 '$sample_id'_trimmed_R1.fastq.gz --out2 '$sample_id'_trimmed_R2.fastq.gz
        """
    else if (params.Library_Layout == "SINGLE")     
        """
        fastp --thread $task.cpus $params.mod_FastP --in1 ${reads[0]} --length_required 30 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 30 --trim_front1=10 --out1 '$sample_id'_trimmed_R1.fastq.gz
        """
    else
        error "Undefined Library_Layout mode: $params.Library_Layout"
}
process FASTP_miRNAs {    //
    maxForks 5    
    publishDir "${params.outdir}/Trimmed", mode:'copy', pattern: '*.fastq.gz' 
    tag "FASTP on $sample_id"    
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads
    
    script:
    if (params.Library_Layout == "PAIRED")    
        """
        fastp --thread $task.cpus $params.mod_FastP --in1 ${reads[0]} --in2 ${reads[1]} --length_required 18 --length_limit 30 --detect_adapter_for_pe --cut_right --cut_right_window_size 4 --cut_right_mean_quality 30 --out1 '$sample_id'_trimmed_R1.fastq.gz --out2 '$sample_id'_trimmed_R2.fastq.gz
        """
    else if (params.Library_Layout == "SINGLE")     
        """
        fastp --thread $task.cpus $params.mod_FastP --in1 ${reads[0]} --length_required 18 --length_limit 30 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 30 --out1 '$sample_id'_trimmed_R1.fastq.gz
        """
    else
        error "Undefined Library_Layout mode: $params.Library_Layout"
}
process strandness_test { 
    conda "/media/storage2/software/anaconda3/envs/MirScience" 
    tag "strandness_test"    
    
    input:
    tuple val(sample_id), path(reads)

    output:
    path "detected_strandedness.txt", emit: strandness_def
        
    script:
    if (params.Library_Layout == "PAIRED")   
        """
        check_strandedness -g $params.GTF --transcripts $params.transcripts_file -r1 ${reads[0]} -r2 ${reads[1]} | tail -1 | cut -d' ' -f 4| cut -d'/' -f 1 > detected_strandedness.txt
        """
    else if (params.Library_Layout == "SINGLE")   
        """
        check_strandedness -g $params.GTF --transcripts $params.transcripts_file -r1 ${reads} | tail -1 | cut -d' ' -f 4| cut -d'/' -f 1 > detected_strandedness.txt
        """
    else
        error "Undefined Library_Layout mode: $params.Library_Layout"
}
process Hisat2 {    //
    maxForks 2
    cpus 40         
    tag "Hisat2 on $sample_id"    
    
    input:
    tuple val(sample_id), path(reads)
    each(strandness)

    output:
    path("${sample_id}.sorted.bam"), emit: bams      
        
    script:
    if (params.Library_Layout == "PAIRED") 
        if (strandness == "RF")
            """
            hisat2 --rna-strandness RF -p $task.cpus --summary-file '$sample_id'.summary -x $params.index_genome -1 ${reads[0]} -2 ${reads[1]} -S '$sample_id'.sam
            java -jar /media/storage2/software/Picard/picard.jar SortSam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate INPUT= '$sample_id'.sam OUTPUT= '$sample_id'.sorted.bam
            """
        else if (strandness == "FR")
            """
            hisat2 --rna-strandness FR -p $task.cpus --summary-file '$sample_id'.summary -x $params.index_genome -1 ${reads[0]} -2 ${reads[1]} -S '$sample_id'.sam
            java -jar /media/storage2/software/Picard/picard.jar SortSam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate INPUT= '$sample_id'.sam OUTPUT= '$sample_id'.sorted.bam
            """
        else if (strandness == "unstranded")
            """
            hisat2 -p $task.cpus --summary-file '$sample_id'.summary -x $params.index_genome -1 ${reads[0]} -2 ${reads[1]} -S '$sample_id'.sam
            java -jar /media/storage2/software/Picard/picard.jar SortSam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate INPUT= '$sample_id'.sam OUTPUT= '$sample_id'.sorted.bam
            """
        else
            error "Undefined strandness mode: $strandness"
    else if (params.Library_Layout == "SINGLE") 
        if (strandness == "RF")
            """
            hisat2 --rna-strandness R -p $task.cpus --summary-file '$sample_id'.summary -x $params.index_genome -U ${reads[0]} -S '$sample_id'.sam
            java -jar /media/storage2/software/Picard/picard.jar SortSam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate INPUT= '$sample_id'.sam OUTPUT= '$sample_id'.sorted.bam
            """
        else if (strandness == "FR")
            """
            hisat2 --rna-strandness F -p $task.cpus --summary-file '$sample_id'.summary -x $params.index_genome -U ${reads[0]} -S '$sample_id'.sam
            java -jar /media/storage2/software/Picard/picard.jar SortSam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate INPUT= '$sample_id'.sam OUTPUT= '$sample_id'.sorted.bam
            """
        else if (strandness == "unstranded")
            """
            hisat2 -p $task.cpus --summary-file '$sample_id'.summary -x $params.index_genome -U ${reads[0]} -S '$sample_id'.sam
            java -jar /media/storage2/software/Picard/picard.jar SortSam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate INPUT= '$sample_id'.sam OUTPUT= '$sample_id'.sorted.bam
            """
        else
            error "Undefined strandness mode: $strandness"
    else
        error "Undefined Library_Layout mode: $params.Library_Layout"
}
process FeatureCounts {    //  
    publishDir params.outdir, mode:'copy'        
    
    input:
    file archivos 
    val(strandness)

    output:
    path("Matrix_count.txt"), emit: matrix_count  
        
    script:
    if (params.Library_Layout == "PAIRED") 
        if (strandness == "RF")
            """
            featureCounts -s 2 -T $task.cpus -p -t exon -g gene_id -a $params.GTF -o matrix_count.txt $archivos
            grep -vE "^#" matrix_count.txt > Matrix_count.txt        
            """
        else if (strandness == "FR")
            """
            featureCounts -s 1 -T $task.cpus -p -t exon -g gene_id -a $params.GTF -o matrix_count.txt $archivos
            grep -vE "^#" matrix_count.txt > Matrix_count.txt
            """
        else if (strandness == "unstranded")
            """
            featureCounts -s 0 -T $task.cpus -p -t exon -g gene_id -a $params.GTF -o matrix_count.txt $archivos
            grep -vE "^#" matrix_count.txt > Matrix_count.txt
            """
        else
            error "Undefined strandness mode: $strandness"
    else if (params.Library_Layout == "SINGLE") 
        if (strandness == "RF")
            """
            featureCounts -s 2 -T $task.cpus -t exon -g gene_id -a $params.GTF -o matrix_count.txt $archivos
            grep -vE "^#" matrix_count.txt > Matrix_count.txt        
            """
        else if (strandness == "FR")
            """
            featureCounts -s 1 -T $task.cpus -t exon -g gene_id -a $params.GTF -o matrix_count.txt $archivos
            grep -vE "^#" matrix_count.txt > Matrix_count.txt
            """
        else if (strandness == "unstranded")
            """
            featureCounts -s 0 -T $task.cpus -t exon -g gene_id -a $params.GTF -o matrix_count.txt $archivos
            grep -vE "^#" matrix_count.txt > Matrix_count.txt
            """
        else
            error "Undefined strandness mode: $strandness"
    else
        error "Undefined Library_Layout mode: $params.Library_Layout"
}
process Process_miRNAs_fastq {
    conda "/media/storage2/software/anaconda3/envs/miRDeep2" 
    maxForks 5       
    tag "Process_miRNAs_fastq on $sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:    
    tuple val(sample_id), path("${sample_id}.tab"), emit: quantified_miRNAs

    script:
    """
    gunzip -c $reads > '$sample_id'.fastq
    fastq2fasta.pl '$sample_id'.fastq > '$sample_id'.fa
    collapse_reads_md.pl '$sample_id'.fa $params.organism_code > '$sample_id'.collapsed.fa
    quantifier.pl -p $params.precursor_miRNAs -m $params.mature_miRNAs -r '$sample_id'.collapsed.fa -t $params.organism_code -d -j
    mv *.csv '$sample_id'.tab
    """
}
process Process_read_counts_files {    
    maxForks 2 

    input:
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}.fixed.tab"), emit: quantified_miRNAs
    
    script:
    """   
    #!/usr/bin/env Rscript 
    expr0 <- read.delim("${reads}", sep = "\t")
    expr0 <- expr0[c("X.miRNA","precursor","read_count")]
    names(expr0)[1] <- "miRNA"
    expr0\$ID <- paste0(expr0\$precursor,"|",expr0\$miRNA)
    expr0 <- expr0[c("ID","read_count")]
    expr0 <- aggregate(x = expr0\$read_count, by = list(expr0\$ID), FUN = sum)
    names(expr0) <- c("ID","read_count")
    write.table(expr0,file="${sample_id}.fixed.tab", quote = F, row.names = F, sep = "\t")
    """
}
process Parse_read_counts_files {
    conda "/media/storage2/software/anaconda3/envs/MirScience"     
    publishDir params.outdir, mode:'copy'  

    input:
    file archivos 

    output:
    path("miRNAs_matrix.tab"), emit: matrix_count  
    
    script:
    """   
    #!/usr/bin/env python3
    import pandas as pd 
    import os
    import sys
    from os import walk
    Dir_files = pd.DataFrame(next(walk("./"), (None, None, []))[2], columns={"file"})
    Sample_files = Dir_files.loc[Dir_files.file.str.contains(".fixed.tab")].copy()

    miRNAs_matrix = pd.DataFrame()
    for sample in Sample_files.file:
        sample_name = sample.split(".")[0]
        sample_df = pd.read_csv(sample,sep="\t",header=0,index_col=0).rename(columns={"read_count":sample_name})
        miRNAs_matrix = pd.merge(miRNAs_matrix,sample_df, left_index=True,right_index=True, how="outer").fillna(0)
        print(sample_df)
    
    miRNAs_matrix.to_csv("miRNAs_matrix.tab", sep = "\t")
    """
}
process DESeq2 {    
    publishDir params.outdir, mode:'copy'  

    input:
    path(Matrix)
    path(Sample_table)
    path(Annotation)
    val(Project)

    output:
    path("*.rds"), emit: rds
    path "*"
    
    script:
    """   
    #!/usr/bin/env Rscript 
    library(pdftools)
    library(DESeq2)
    library(ggplot2)
    library(BiocParallel)
    library(stringr)
    register(MulticoreParam(15))
    expr0 <- read.delim("${Matrix}", sep = "\t")
    colnames(expr0) <- gsub(".", "-", colnames(expr0),fixed=T)
    colnames(expr0) <- gsub("-sorted-bam", "", colnames(expr0))
    expr0 <- expr0[!grepl("_PAR_Y",expr0\$Geneid),]
    rownames(expr0) <- sapply(strsplit(expr0\$Geneid,".",fixed=T), "[", 1)
    expr0\$Geneid <- NULL
    expr0\$Chr <- NULL
    expr0\$Length <- NULL
    expr0\$Start <- NULL
    expr0\$End <- NULL
    expr0\$Strand <- NULL
    gencode_info <- read.delim("${Annotation}", sep = "\t")
    gencode_info\$gene_id <- sapply(strsplit(gencode_info\$gene_id,".",fixed=T), "[", 1)
    sample_data <- read.delim("${Sample_table}", sep = ",")
    sample_data <- sample_data[sample_data\$BioProjectID == "${Project}",]
    sample_data <- sample_data[sample_data\$Condition != "ND",]    
    sample_data\$Condition <- gsub("1", "${params.Case_name}", sample_data\$Condition)
    sample_data\$Condition <- gsub("0", "${params.Control_name}", sample_data\$Condition) 
    write.csv(sample_data, file="Samples_${Project}_${params.Control_name}_vs_${params.Case_name}.csv", row.names = F)
    sample_data <- sample_data[c("Run_acc","Condition")]
    names(sample_data) <- c("column", "Condition")
    coldata <- sample_data 
    matrix <- as.data.frame(expr0)    
    RNAs <- "${Project}"   
    cts <- matrix[coldata\$column]     
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)
    keep <- rowSums(counts(dds)) > 1 ## Pre-Filtering
    dds <- dds[keep,]
    dds\$Condition <- factor(dds\$Condition, levels = c("${params.Control_name}","${params.Case_name}"))
    dds\$Condition <- droplevels(dds\$Condition)
    dds <- DESeq(dds, fitType="local", parallel = TRUE)        
    DESeq_norm <- counts(dds, normalized=T)
    write.table(DESeq_norm, file = paste0(RNAs,"_${Project}_${params.Control_name}_vs_${params.Case_name}_DESeq_norm_matrix_WT.csv"), row.names = T, quote = F, col.names = T, sep = ",")  
    vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    pcaData <- plotPCA(vsd, intgroup="Condition", returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    p <- ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
    ggtitle(paste(RNAs, sep = "")) +          
    theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold")) +
    geom_point(size=2) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    scale_color_manual(values=c("#56B4E9", "red"))
    ggsave(paste0(RNAs,"_${Project}_${params.Control_name}_vs_${params.Case_name}_PCA_plot.png"), plot = p, width = 8, height = 8, dpi = 300, units = "in")
    res <- results(dds, parallel = TRUE)
    res2 <- as.data.frame(res)
    res2 <- merge(res2, gencode_info, by.x = "row.names",by.y = "gene_id")
    res2\$L2FC_ABS <- abs(res2\$log2FoldChange)
    res2 <- res2[order(res2\$L2FC_ABS,decreasing=T),]
    write.csv(na.omit(res2[res2\$padj < 0.05,]), file=paste(RNAs,"_",levels(dds\$Condition)[1], "_vs_", levels(dds\$Condition)[2],"_DE.csv", sep= ""))
    write.csv(res2, file=paste(RNAs,"_",levels(dds\$Condition)[1], "_vs_", levels(dds\$Condition)[2],"_DE_NS.csv", sep= ""))
    saveRDS(na.omit(res2[res2\$padj < 0.05,]), file = "${Project}_DEGs.rds")
    """
}
process DESeq2_miRNAs {    
    publishDir params.outdir, mode:'copy'  

    input:
    path(Matrix)
    path(Sample_table)
    val(Project)

    output:
    path "*"
    
    script:
    """   
    #!/usr/bin/env Rscript 
    library(pdftools)
    library(DESeq2)
    library(ggplot2)
    library(BiocParallel)
    library(stringr)
    register(MulticoreParam(15))
    expr0 <- read.delim("${Matrix}", sep = "\t")
    colnames(expr0) <- gsub(".", "-", colnames(expr0),fixed=T)
    rownames(expr0) <-expr0\$ID
    sample_data <- read.delim("${Sample_table}", sep = ",")
    sample_data <- sample_data[sample_data\$BioProjectID == "${Project}",]
    sample_data <- sample_data[sample_data\$Condition != "ND",]
    sample_data\$Condition <- gsub("1", "${params.Case_name}", sample_data\$Condition)
    sample_data\$Condition <- gsub("0", "${params.Control_name}", sample_data\$Condition) 
    write.csv(sample_data, file="Samples_${Project}_${params.Control_name}_vs_${params.Case_name}.csv", row.names = F)
    sample_data <- sample_data[c("Run_acc","Condition")]
    names(sample_data) <- c("column", "Condition")
    coldata <- sample_data 
    matrix <- as.data.frame(expr0)    
    RNAs <- "${Project}"   
    cts <- matrix[coldata\$column]     
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)
    keep <- rowSums(counts(dds)) > 1 ## Pre-Filtering
    dds <- dds[keep,]
    dds\$Condition <- factor(dds\$Condition, levels = c("${params.Control_name}","${params.Case_name}"))
    dds\$Condition <- droplevels(dds\$Condition)
    dds <- DESeq(dds, fitType="local", parallel = TRUE)        
    DESeq_norm <- counts(dds, normalized=T)
    write.table(DESeq_norm, file = paste0(RNAs,"_${Project}_${params.Control_name}_vs_${params.Case_name}_DESeq_norm_matrix_WT.csv"), row.names = T, quote = F, col.names = T, sep = ",")  
    vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    pcaData <- plotPCA(vsd, intgroup="Condition", returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    p <- ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
    ggtitle(paste(RNAs, sep = "")) +          
    theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold")) +
    geom_point(size=2) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    scale_color_manual(values=c("#56B4E9", "red"))
    ggsave(paste0(RNAs,"_${Project}_${params.Control_name}_vs_${params.Case_name}_PCA_plot.png"), plot = p, width = 8, height = 8, dpi = 300, units = "in")
    res <- results(dds, parallel = TRUE)
    res2 <- as.data.frame(res)
    res2\$L2FC_ABS <- abs(res2\$log2FoldChange)
    res2 <- res2[order(res2\$L2FC_ABS,decreasing=T),]
    write.csv(na.omit(res2[res2\$padj < 0.05,]), file=paste(RNAs,"_",levels(dds\$Condition)[1], "_vs_", levels(dds\$Condition)[2],"_DE.csv", sep= ""))
    write.csv(res2, file=paste(RNAs,"_",levels(dds\$Condition)[1], "_vs_", levels(dds\$Condition)[2],"_DE_NS.csv", sep= ""))
    """
}
process Enrichr { 
    publishDir "${params.outdir}", mode:'copy'
    maxForks 1   
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(rds_file)

    output:  
    path("enrichment.pdf")

    script:
    """
    #!/usr/bin/env Rscript 
    library(enrichR)
    important_genes <- readRDS("${rds_file}")
    genes_alterados <- unique(important_genes\$gene_name)
    setEnrichrSite("Enrichr") 
    websiteLive <- TRUE
    dbs <- listEnrichrDbs()
    dbs <- c("$params.enrich_gene_set")
    enriched <- enrichr(genes_alterados, dbs)   
    enriched <- enriched[[1]][enriched[[1]]\$Adjusted.P.value < 0.05,]
    pdf("enrichment.pdf",width = 16, height=9)
    plotEnrich(enriched, showTerms = 10, numChar = 40, y = "Count", orderBy = "P.value", title = paste(dbs[1], " enrichment\n", sep = "" ))
    dev.off()
    """   
}
workflow QC1{
    take: data
    main:
        FASTQC(data) 
        MULTIQC(FASTQC.out.multiqc_input.collect(),"QC1")
}
workflow Trimm_QC2{
    take: data
    main:
        FASTP(data)    
        FASTQC(FASTP.out.trimmed_reads)    
        MULTIQC(FASTQC.out.multiqc_input.collect(),"QC2")
        emit:
            FASTP.out.trimmed_reads
}
workflow Trimm_miRNAs{
    take: data
    main:
        FASTP_miRNAs(data)    
        FASTQC(FASTP_miRNAs.out.trimmed_reads)    
        MULTIQC(FASTQC.out.multiqc_input.collect(),"QC2")
        emit:
            FASTP_miRNAs.out.trimmed_reads
}
workflow Trimm_short_reads{
    take: data
    main:
        FASTP_short_reads(data)    
        FASTQC(FASTP_short_reads.out.trimmed_reads)    
        MULTIQC(FASTQC.out.multiqc_input.collect(),"QC2")
        emit:
            FASTP_short_reads.out.trimmed_reads
}
workflow Map_and_count {
    take: data
    main:   
        //Trimmed_ch = channel.fromFilePairs(data, checkIfExists: false)
        strandness_test(data.take(1)) 
        strandness_def = strandness_test.out.strandness_def.splitText{it.strip()}
        strandness_def.view()
        Hisat2(data,strandness_def)
        FeatureCounts(Hisat2.out.bams.collect(),strandness_def)
    emit: 
        FeatureCounts.out.matrix_count 
}
workflow Map_and_count_miRNAs {
    take: data
    main:   
        Process_miRNAs_fastq(data)
        Process_read_counts_files(Process_miRNAs_fastq.out.quantified_miRNAs)
        Parse_read_counts_files(Process_read_counts_files.out.quantified_miRNAs.collect())
    emit: 
        Parse_read_counts_files.out.matrix_count 
}
workflow { 
    params.outdir = "Resultados/${params.Organism}/${params.Library_Layout}/${params.Project}/"
    out_dir = file(params.outdir)
    out_dir.mkdir()
    ST = Channel.fromPath(params.sample_table).splitCsv(header: true).filter(row -> row.Condition != "ND").filter(row -> row.BioProjectID == params.Project)
    Sample_table = Channel.fromPath(params.sample_table)
    Projects = ST.map { row-> tuple(row.Organism, row.Library_Layout, row.BioProjectID, row.Platform, row.Run_acc) } 
    Info = Projects.map { it[0,1,2,3] }.unique()
    Runs = Projects.map { it[4] }
    if (params.Organism == "Mus_musculus"){
        params.GTF = "/media/storage/datasets/gtf/mouse/gencode.vM27.annotation.gtf"
        params.Annotation = "/media/storage/Adolfo/MirScience/IDs_equivalencias_gencode.vM27.txt"
        params.transcripts_file = "/media/storage/Adolfo/MirScience/gencode.vM27.transcripts.fa"
        params.index_genome = "/media/storage/datasets/indexes/hisat2_index/Mus_musculus/Mm_GRCm39/GRCm39"
        params.enrich_gene_set= "WikiPathways_2019_Mouse"
    } else if (params.Organism == "Homo_sapiens"){
        params.GTF = "/media/storage/datasets/gtf/human/gencode.v38.annotation.gtf"
        params.Annotation = "/media/storage/Adolfo/MirScience/IDs_equivalencias_gencode.v38.txt"
        params.transcripts_file = "/media/storage/Adolfo/MirScience/gencode.v38.transcripts.fa"
        params.index_genome = "/media/storage/datasets/indexes/hisat2_index/Homo_sapiens/grch38/genome"
        params.enrich_gene_set= "GO_Biological_Process_2021"
    } else {
        close
    }
    if (params.mod_FastP == ""){    
        print "\nCheck Information:\n $params.Project, $params.Library_Layout, NOT Nova or Nextseq, $params.Organism\n"
    } else if (params.mod_FastP == "-g"){
         print "\nCheck Information:\n$params.Project, $params.Library_Layout, IS Nova or Nextseq, $params.Organism\n"}
    print 'If something is wrong press "ctrl + C" and change parameters eg. --mod_FastP="-g" (for novaseq or nextseq data)\nothers params are: Project, Library_layout, Organism'
    print "Workflow project outdir " + params.outdir
    Dowload_fastqs(Runs)
    QC1(Dowload_fastqs.out.fastqs)
    Trimm_QC2(Dowload_fastqs.out.fastqs)
    Map_and_count(Trimm_QC2.out)  
    Annotation_ch = Channel.fromPath(params.Annotation)
    DESeq2(Map_and_count.out,Sample_table,Annotation_ch,params.Project) 
    Enrichr(DESeq2.out.rds) 
}
workflow miRNAseq {
    params.outdir = "Resultados/miRNAseq/${params.Organism}/${params.Library_Layout}/${params.Project}/"
    out_dir = file(params.outdir)
    out_dir.mkdir()
    ST = Channel.fromPath(params.sample_table).splitCsv(header: true).filter(row -> row.Condition != "ND").filter(row -> row.BioProjectID == params.Project)
    Sample_table = Channel.fromPath(params.sample_table)
    Projects = ST.map { row-> tuple(row.Organism, row.Library_Layout, row.BioProjectID, row.Platform, row.Run_acc) } 
    Info = Projects.map { it[0,1,2,3] }.unique()
    Runs = Projects.map { it[4] }
    if (params.Organism == "Mus_musculus"){
        params.mature_miRNAs = "/media/storage/Adolfo/MirScience/RNA-seq/mature_mmu_ref.fa"
        params.precursor_miRNAs = "/media/storage/Adolfo/MirScience/RNA-seq/hairpin_mmu_ref.fa"
        params.organism_code = "mmu"
    } else if (params.Organism == "Homo_sapiens"){
        params.mature_miRNAs = "/media/storage/Adolfo/MirScience/RNA-seq/mature_hsa_ref.fa"
        params.precursor_miRNAs = "/media/storage/Adolfo/MirScience/RNA-seq/hairpin_hsa_ref.fa"
        params.organism_code = "hsa"
    } else {
        close
    }
    if (params.mod_FastP == ""){    
        print "\nCheck Information:\n $params.Project, $params.Library_Layout, NOT Nova or Nextseq, $params.Organism\n"
    } else if (params.mod_FastP == "-g"){
         print "\nCheck Information:\n$params.Project, $params.Library_Layout, IS Nova or Nextseq, $params.Organism\n"}
    print 'If something is wrong press "ctrl + C" and change parameters eg. --mod_FastP="-g" (for novaseq or nextseq data)\nothers params are: Project, Library_layout, Organism'
    print "Workflow project outdir " + params.outdir
    Dowload_fastqs(Runs)
    QC1(Dowload_fastqs.out.fastqs)
    Trimm_miRNAs(Dowload_fastqs.out.fastqs)
    Map_and_count_miRNAs(Trimm_miRNAs.out)
    DESeq2_miRNAs(Map_and_count_miRNAs.out,Sample_table,params.Project)  
}
workflow short_reads_RNAseq { 
    params.outdir = "Resultados/${params.Organism}/${params.Library_Layout}/${params.Project}/"
    out_dir = file(params.outdir)
    out_dir.mkdir()
    ST = Channel.fromPath(params.sample_table).splitCsv(header: true).filter(row -> row.Condition != "ND").filter(row -> row.BioProjectID == params.Project)
    Sample_table = Channel.fromPath(params.sample_table)
    Projects = ST.map { row-> tuple(row.Organism, row.Library_Layout, row.BioProjectID, row.Platform, row.Run_acc) } 
    Info = Projects.map { it[0,1,2,3] }.unique()
    Runs = Projects.map { it[4] }
    if (params.Organism == "Mus_musculus"){
        params.GTF = "/media/storage/datasets/gtf/mouse/gencode.vM27.annotation.gtf"
        params.Annotation = "/media/storage/Adolfo/MirScience/IDs_equivalencias_gencode.vM27.txt"
        params.transcripts_file = "/media/storage/Adolfo/MirScience/gencode.vM27.transcripts.fa"
        params.index_genome = "/media/storage/datasets/indexes/hisat2_index/Mus_musculus/Mm_GRCm39/GRCm39"
        params.enrich_gene_set= "WikiPathways_2019_Mouse"
    } else if (params.Organism == "Homo_sapiens"){
        params.GTF = "/media/storage/datasets/gtf/human/gencode.v38.annotation.gtf"
        params.Annotation = "/media/storage/Adolfo/MirScience/IDs_equivalencias_gencode.v38.txt"
        params.transcripts_file = "/media/storage/Adolfo/MirScience/gencode.v38.transcripts.fa"
        params.index_genome = "/media/storage/datasets/indexes/hisat2_index/Homo_sapiens/grch38/genome"
        params.enrich_gene_set= "GO_Biological_Process_2021"
    } else {
        close
    }
    if (params.mod_FastP == ""){    
        print "\nCheck Information:\n $params.Project, $params.Library_Layout, NOT Nova or Nextseq, $params.Organism\n"
    } else if (params.mod_FastP == "-g"){
         print "\nCheck Information:\n$params.Project, $params.Library_Layout, IS Nova or Nextseq, $params.Organism\n"}
    print 'If something is wrong press "ctrl + C" and change parameters eg. --mod_FastP="-g" (for novaseq or nextseq data)\nothers params are: Project, Library_layout, Organism'
    print "Workflow project outdir " + params.outdir
    Dowload_fastqs(Runs)
    QC1(Dowload_fastqs.out.fastqs)
    Trimm_short_reads(Dowload_fastqs.out.fastqs)
    Map_and_count(Trimm_short_reads.out)  
    Annotation_ch = Channel.fromPath(params.Annotation)
    DESeq2(Map_and_count.out,Sample_table,Annotation_ch,params.Project) 
    Enrichr(DESeq2.out.rds) 
}
workflow Download_only { 
    params.outdir = "Resultados/${params.Organism}/${params.Library_Layout}/${params.Project}/"
    out_dir = file(params.outdir)
    out_dir.mkdir()
    ST = Channel.fromPath(params.sample_table).splitCsv(header: true).filter(row -> row.Condition != "ND").filter(row -> row.BioProjectID == params.Project)
    Sample_table = Channel.fromPath(params.sample_table)
    Projects = ST.map { row-> tuple(row.Organism, row.Library_Layout, row.BioProjectID, row.Platform, row.Run_acc) } 
    Info = Projects.map { it[0,1,2,3] }.unique()
    Runs = Projects.map { it[4] }
    if (params.mod_FastP == ""){    
        print "\nCheck Information:\n $params.Project, $params.Library_Layout, NOT Nova or Nextseq, $params.Organism\n"
    } else if (params.mod_FastP == "-g"){
         print "\nCheck Information:\n$params.Project, $params.Library_Layout, IS Nova or Nextseq, $params.Organism\n"}
    print 'If something is wrong press "ctrl + C" and change parameters eg. --mod_FastP="-g" (for novaseq or nextseq data)\nothers params are: Project, Library_layout, Organism'
    print "Workflow project outdir " + params.outdir
    Dowload_fastqs(Runs)
}
workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
    //def msg = """\
    //    Pipeline execution summary
    //    ---------------------------
    //    Completed at: ${workflow.complete}
    //    Duration    : ${workflow.duration}
    //    Success     : ${workflow.success}
    //    workDir     : ${workflow.workDir}
    //    exit status : ${workflow.exitStatus}
    //    """
    //    .stripIndent()

    //sendMail(to: 'adolfo.rojas@ug.uchile.cl', subject: 'My pipeline execution', body: msg)
}
