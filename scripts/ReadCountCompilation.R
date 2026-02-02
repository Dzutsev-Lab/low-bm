fastq_counts <- read.delim(file = stdin(), header =  TRUE)
dada_counts <- read.delim(snakemake@input$dada_counts, header =TRUE)
seq_table <- read.delim(snakemake@input$seq_table, header = TRUE, row.names = 1)
host_unmapped_ASV_names <- readLines(snakemake@input$unmapped_host_ASV_names)
viral_unmapped_ASV_names <- readLines(snakemake@input$unmapped_viral_ASV_names)
bacterial_mapped_ASV_names <- readLines(snakemake@input$bacterial_ASV_names)

host_unmapped_seq_table <- seq_table[, host_unmapped_ASV_names]
viral_unmapped_seq_table <- seq_table[, viral_unmapped_ASV_names]
bacterial_mapped_seq_table <- seq_table[, bacterial_mapped_ASV_names]



host_unmapped_seq_table$HostUnmapped <- rowSums(host_unmapped_seq_table)
host_unmapped_seq_table$SampleID <- rownames(host_unmapped_seq_table)

viral_unmapped_seq_table$ViralUnmapped <- rowSums(viral_unmapped_seq_table)
viral_unmapped_seq_table$SampleID <- rownames(viral_unmapped_seq_table)

bacterial_mapped_seq_table$BacterialMapped <- rowSums(bacterial_mapped_seq_table)
bacterial_mapped_seq_table$SampleID <- rownames(bacterial_mapped_seq_table)

combined_counts <- full_join(fastq_counts, dada_counts, host_unmapped_seq_table, viral_unmapped_seq_table, bacterial_mapped_seq_table, by = "SampleID")


