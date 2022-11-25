setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(GenomicRanges)
require(readxl)

#----------------------------------------------------------#
# Import metadata
#----------------------------------------------------------#
meta <- as.data.table(read_xlsx("Rdata/metadata_RNA.xlsx"))
# meta <- meta[!is.na(DESeq2_object)]
meta <- meta[DESeq2_object %in% c("epiCancer_GFP", "epiCancer_noGFP")]
# Retrieve fq files 
meta[, c("fq1", "fq2"):= lapply(.SD, function(x){
  list.files(paste0("/mnt/f/_R_data/projects/epigenetic_cancer/db/fastq/", project), 
             paste0("^", x, "$"), 
             recursive= T,
             full.names = T)
}), by= .(fq1, fq2, project), .SDcols= c("fq1", "fq2")]

#----------------------------------------------------------#
# BUILD dm6 index with GFP sequences
#----------------------------------------------------------#
# Build dm6 index
if(!file.exists("/mnt/d/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index.log"))
{
  ref <- seqinr::read.fasta("/mnt/d/_R_data/genomes/dm6/Sequence/WholeGenomeFasta/genome.fa",
                            as.string = T)
  ref <- lapply(ref, as.character)
  ref <- c(ref, 
           GFP= "gtgagcaagggcgagaagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagatgtccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcaaaaccaccctgacctggggcatgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcgtcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaacgccatcagcggcaacgccaatatcaccgccgacaagcagaagaacggcatcaaggcctacttcacgatccgccacgacgtcgaggacggcagcgtgctgctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccaagcagagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatccctctcggcgcggacgagctgtacaag",
           EGFP= "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA",
           mRFP1= "ATGGCCTCCTCCGAGGACGTCATCAAGGAGTTCATGCGCTTCAAGGTGCGCATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGCGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCCAGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCACCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGATGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCCGAGGTCAAGACCACCTACATGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAAGACCGACATCAAGCTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAGCGCGCCGAGGGCCGCCACTCCACCGGCGCCTAA")
  ref <- lapply(ref, tolower)
  seqinr::write.fasta(sequences = ref, 
                      names= names(ref), 
                      "/mnt/d/_R_data/genomes/dm6_GFP/dm6_GFP.fa", 
                      as.string = T)
  buildindex(basename= "/mnt/d/_R_data/genomes/dm6_GFP/subreadr_index/subreadr_dm6_GFP_index", 
             reference= "/mnt/d/_R_data/genomes/dm6_GFP/dm6_GFP.fa")
}

#----------------------------------------------------------#
# ALIGNMENT
#----------------------------------------------------------#
meta[, bam:= paste0("/mnt/f/_R_data/projects/epigenetic_cancer/db/bam/", project, "/",
                    DESeq2_object, "_", cdition, "_rep", rep, ".bam"), .(DESeq2_object, cdition, rep)]
meta[, {
  if(!file.exists(bam))
  {
    print(paste0("START...", bam))
    align(index= "/mnt/d/_R_data/genomes/dm6_GFP/subreadr_index/subreadr_dm6_GFP_index",
          readfile1= fq1,
          readfile2= ifelse(is.na(fq2), NULL, fq2),
          maxMismatches = 6,
          nthreads = 10,
          unique = T,
          output_file= bam)
  }
  print(paste(bam, "--> DONE!"))
}, .(bam, fq1, fq2)]

#----------------------------------------------------------#
# GFP SAF file
#----------------------------------------------------------#
fwrite(data.table(GeneID= c("GFP", "EGFP", "mRFP1"), 
                  Chr=  c("GFP", "EGFP", "mRFP1"),
                  start= 1,
                  end= 800,
                  Strand= "+"),
       file = "db/saf/GFP_genes.saf",
       col.names = T,
       sep= "\t",
       quote= F)

#----------------------------------------------------------#
# Counts
#----------------------------------------------------------#
meta[, read_counts:= paste0("db/counts/RNA/", DESeq2_object, "_read_counts.rds"), DESeq2_object]
meta[, GFP_counts:= paste0("db/counts/RNA/", DESeq2_object, "_GFP_counts.rds"), DESeq2_object]
meta[, {
  if(!file.exists(read_counts))
  {
    .c <- featureCounts(bam, # count reads
                        annot.ext= "../../genomes/dm6/dmel-all-r6.36.gtf",
                        isGTFAnnotationFile = T,
                        isPairedEnd = layout=="PAIRED",
                        nthreads = 8)
    saveRDS(.c, read_counts)
  }
  if(!file.exists(GFP_counts))
  {
    .c <- featureCounts(bam, # count reads
                        annot.ext= "db/saf/GFP_genes.saf",
                        isGTFAnnotationFile = F,
                        isPairedEnd = layout=="PAIRED",
                        nthreads = 8)
    saveRDS(.c, GFP_counts)
  }
}, .(read_counts, GFP_counts, layout)]

#----------------------------------------------------------#
# DESeq2
#----------------------------------------------------------#
meta[, dds_file:= paste0("db/dds/RNA/", DESeq2_object, "_dds.rds"), DESeq2_object]
FC <- meta[, {
  if(!file.exists(dds_file))
  {
    # SampleTable and DF
    .c <- readRDS(read_counts)
    DF <- as.data.frame(.c[[1]])
    DF <- DF[order(rownames(DF)),]
    colnames(DF) <- paste0(cdition, "_rep", rep)
    sampleTable <- as.data.frame(setNames(tstrsplit(names(DF), "_"),
                                          c("cdition", "rep")),
                                 row.names = names(DF))
    # Filter low read genes
    DF <- DF[rowSums(DF)>10, ] # filter low read genes + N
    # dds
    dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                          colData= sampleTable,
                                          design= ~rep+cdition)
    dds <- DESeq2::DESeq(dds)
    # Add gene length to compute fpkms
    glength <- as.data.table(.c$annotation, key= "GeneID")[, .(GeneID, Length)]
    mcols(dds)$basepairs <- glength[names(dds@rowRanges), Length]
    # SAVE
    saveRDS(dds, dds_file)
  }else
    dds <- readRDS(dds_file)
  CJ(as.character(dds$cdition),
     as.character(dds$cdition), 
     unique = T)[V1!=V2]
}, .(read_counts, dds_file, DESeq2_object)]

#----------------------------------------------------------#
# FC tables
#----------------------------------------------------------#
cditions_table <- read_xlsx("Rdata/RNA_dds_conditions.xlsx")
setkeyv(FC, c("V1", "V2"))
FC <- FC[as.data.table(cditions_table), on= c("V1==condition", "V2==control"), nomatch= NULL]
FC[, FC_file:= paste0("db/FC_tables/RNA/", DESeq2_object, "_", V1, "_vs_", V2, ".txt"), .(DESeq2_object, V1, V2)]
FC[, {
  dds <- readRDS(dds_file)
  .SD[, {
    if(!file.exists(FC_file))
    {
      # Compute FC tables
      .c <- DESeq2::lfcShrink(dds,
                              type= "ashr",
                              contrast= c("cdition", V1, V2))
      .c <- as.data.table(as.data.frame(.c), 
                          keep.rownames = "FBgn")
      .c[, diff:= fcase(padj<0.05 & log2FoldChange>1, "up",
                        padj<0.05 & log2FoldChange<(-1), "down",
                        default= "unaffected")]
      fwrite(.c,
             FC_file, 
             col.names = T,
             sep= "\t")
    }
    print("DONE")
  }, .(V1, V2, FC_file)]
  print("DONE")
}, dds_file]
# Add to meta
meta[FC, FC_file:= i.FC_file, on= c("cdition==V1", "dds_file")]

#----------------------------------------------------------#
# BW files
#----------------------------------------------------------#
meta[DESeq2_object=="epiCancer_noGFP", bw_merge:= paste0("db/bw/RNA/", DESeq2_object, "_", cdition, "_merge.bw")]
meta[DESeq2_object=="epiCancer_noGFP", {
  if(!file.exists(bw_merge))
  {
    bed <- fread(cmd= paste0("bamToBed -i ", bam),
                 fill= T,
                 sep= "\t",
                 sel= c(1,2,3,5),
                 col.names= c("seqnames", "start", "end", "mapq"))
    bed <- GRanges(bed)
    cov <- coverage(bed)/length(bed)*1e6
    rtracklayer::export.bw(GRanges(cov),
                           con= bw_merge)
  }
  print("done")
}, bw_merge]

#----------------------------------------------------------#
# SAVE
#----------------------------------------------------------#
meta[DESeq2_object=="epiCancer_GFP", system:= "GFP"]
meta[DESeq2_object=="epiCancer_noGFP", system:= "noGFP"]
fwrite(meta, 
       "Rdata/processed_metadata_RNA.txt",
       na = NA)
