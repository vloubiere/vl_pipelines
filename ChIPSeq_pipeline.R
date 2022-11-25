setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)
require(data.table)
require(parallel)
require(Rsubread)
require(seqinr)
require(rtracklayer)

#--------------------------------------------------------------#
# METDATA
#--------------------------------------------------------------#
meta <- readxl::read_xlsx("Rdata/metadata_cutnrun.xlsx")
meta <- as.data.table(meta)[!Comment %in% c("failed", "test")]
meta[is.na(Suffix), Suffix:= ""]

#--------------------------------------------------------------#
# Scer/Dm6 combined index
# Yeast genome -> https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna
#--------------------------------------------------------------#
if(length(list.files("/mnt/d/_R_data/genomes/dm6_S288C_combined_bowtie2/", ".bt2$"))==0)
{
  dm6 <- read.fasta("/mnt/d/_R_data/genomes/dm6/Sequence/WholeGenomeFasta/genome.fa")
  S288C <- read.fasta("/mnt/d/_R_data/genomes/S288C/GCF_000146045.2_R64_genomic.fna")
  cmb <- c(dm6, S288C)
  write.fasta(sequences = cmb, 
              names = names(cmb),
              file.out = "/mnt/d/_R_data/genomes/dm6_S288C_combined/dm6_S288C_combined.fa")
  system("bowtie2-build /mnt/d/_R_data/genomes/dm6_S288C_combined_bowtie2/dm6_S288C.fna /mnt/d/_R_data/genomes/dm6_S288C_combined_bowtie2/dm6_S288C")
}

#--------------------------------------------------------------#
# Download fastq files using sra toolkit
#--------------------------------------------------------------#
meta[, vl_bsub(paste("module load build-env/.f2021; module load build-env/.f2021; module load sra-toolkit/3.0.0-centos_linux64;
                     fastq-dump -A", Run, "--split-files --gzip -O /groups/stark/vloubiere/projects/insulator_tomas/db/fastq/")), Run]

#--------------------------------------------------------------#
# Retrieve fastqs and trim reads
#--------------------------------------------------------------#
meta[, fq1:= list.files("/mnt/f/_R_data/projects/epigenetic_cancer/db/fastq/Cut_n_run/", 
                        recursive = T, 
                        full.names = T, 
                        pattern = paste0(fq1, "$")), fq1]
meta[, fq2:= list.files("/mnt/f/_R_data/projects/epigenetic_cancer/db/fastq/Cut_n_run/", 
                        recursive = T, 
                        full.names = T, 
                        pattern = paste0(fq2, "$")), fq2]

meta[, fq1_trim:= gsub(".fq.gz$", "_val_1.fq.gz", fq1)]
meta[, fq2_trim:= gsub(".fq.gz$", "_val_2.fq.gz", fq2)]

mcmapply(FUN = function(outdir, fq1, fq2, fq1_trim, fq2_trim){
  if(!file.exists(fq1_trim) && is.na(fq2_trim))
  {
    cmd <- paste0("trim_galore --paired --gzip -o ", outdir, "/ ", fq1)
    cmd <- paste(c("module load build-env/2020", "module load trim_galore/0.6.2-foss-2018b-python-3.6.6", cmd), collapse= ";")
    system(cmd)
  }
  if((!file.exists(fq1_trim) | !file.exists(fq2_trim)) && !is.na(fq2_trim))
  {
    cmd <- paste0("trim_galore --paired --gzip -o ", outdir, "/ ", fq1, " ", fq2)
    cmd <- paste(c("module load build-env/2020", "module load trim_galore/0.6.2-foss-2018b-python-3.6.6", cmd), collapse= ";")
    system(cmd)
  }
  print("done")
},
outdir= dirname(meta$fq1),
fq1= meta$fq1,
fq2= meta$fq2,
fq1_trim= meta$fq1_trim,
fq2_trim= meta$fq2_trim,
mc.preschedule = F,
mc.cores = getDTthreads())

#--------------------------------------------------------------#
# Alignment
#--------------------------------------------------------------#
# Make chrom_sizes object
chrom_sizes <- rbind(data.table(fread("/mnt/d/_R_data/genomes/dm6/dm6.chrom.sizes.txt", 
                                      col.names = c("seqnames", "seqlengths")), 
                                type="ChIP"),
                     data.table(fread("/mnt/d/_R_data/genomes/S288C/S288C_contigs.txt"),
                                type="spike"))
setkeyv(chrom_sizes, "type")
meta[, bam:= paste0("/mnt/f/_R_data/projects/epigenetic_cancer/db/bam/cutnrun/", ChIP, "_", cdition, "_", rep, Suffix, ".bam")]
meta[, {
  if(!file.exists(bam))
  {
    # bowtie 2
    sam_file <- gsub(".bam$", ".sam", bam)
    cmd <- paste0("bowtie2 -p 10 -x /mnt/d/_R_data/genomes/dm6_S288C_combined_bowtie2/dm6_S288C --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700")
    cmd <- paste0(cmd, " -1 ", fq1_trim)
    cmd <- paste0(cmd, " -2 ", fq2_trim)
    cmd <- paste0(cmd, " -S ", sam_file)
    system(cmd)
    # sam to bam
    cmd <- paste0("/usr/bin/samtools view -@ 9 -b -q 30 ", sam_file, " -o ", bam)
    system(cmd)
    file.remove(sam_file)
  }
  print("DONE")
}, bam]

#--------------------------------------------------------------#
# Peak calling
#--------------------------------------------------------------#
meta[!ChIP %in% c("IgG", "input"), input:= ifelse(ChIP=="PH", "input", "IgG")]
meta[!is.na(input), input_bam:= meta[.BY, bam, on= c("ChIP==input", "rep", "cdition")], .(input, rep, cdition)]
MACS <- function(bam_ChIP,
                 bam_Input,
                 output,
                 broad,
                 SPMR)
{
  cmd <- paste0("/home/vloubiere/.local/bin/macs2 callpeak -g dm --keep-dup 1 -f BAMPE",
                ifelse(SPMR, " -B --SPMR", ""),
                " --outdir ", paste0(dirname(output), "/"),
                " -t ", paste(bam_ChIP, collapse= " "), 
                " -c ", paste(bam_Input, collapse= " "), 
                " -n ", basename(output))
  if(broad)
    cmd <- paste0(cmd, " --broad")
  return(cmd)
}
# Separated replicates
meta[!is.na(input), peaks_rep:= paste0('db/peaks/cutnrun/', paste0(ChIP, "_", cdition, "_", rep))]
meta[!is.na(input), broad:= ChIP %in% c("H3K27me3", "H2AK118Ub", "H3K4me1")]
meta[!is.na(input), {
  if(!file.exists(paste0(peaks_rep, "_peaks.xls")))
  {
    cmd <- MACS(bam, 
                input_bam, 
                peaks_rep, 
                broad = broad,
                SPMR= T)
    system(cmd)
    print("DONE")
  }
}, .(ChIP, peaks_rep, broad)]
# Merge
meta[!is.na(input), peaks_merge:= paste0('db/peaks/cutnrun/', paste0(ChIP, "_", cdition, "_merge"))]
meta[!is.na(input), {
  if(!file.exists(paste0(peaks_merge, "_peaks.xls")))
  {
    cmd <- MACS(bam, 
                input_bam, 
                peaks_merge, 
                broad = broad,
                SPMR= F)
    system(cmd)
    print("DONE")
  }
}, .(ChIP, peaks_merge, broad)]
# Return peaks file
meta[!is.na(input), peaks_rep:= paste0(peaks_rep, ifelse(broad, "_peaks.broadPeak", "_peaks.narrowPeak"))]
meta[!is.na(input), peaks_merge:= paste0(peaks_merge, ifelse(broad, "_peaks.broadPeak", "_peaks.narrowPeak"))]

#--------------------------------------------------------------#
# Confident peaks
#--------------------------------------------------------------#
meta[!is.na(input), filtered_peaks:= paste0("db/peaks/cutnrun/", ChIP, "_", cdition, 
                                            ifelse(broad, 
                                                   "_confident_peaks.broadPeak", 
                                                   "_confident_peaks.narrowPeak"))]
meta[ChIP=="H3K27me3", dist_cutoff:= 2500]
meta[ChIP=="H2AK118Ub", dist_cutoff:= 2500]
meta[ChIP=="H3K27Ac", dist_cutoff:= 250]
meta[ChIP=="H3K36me3", dist_cutoff:= 250]
meta[ChIP=="H3K4me1", dist_cutoff:= 250]
meta[ChIP=="PH", dist_cutoff:= 250]
meta[!is.na(input), {
  if(!file.exists(filtered_peaks))
  {
    .c <- vl_importBed(peaks_merge)
    .c <- .c[signalValue>2 & qValue>2]
    .c <- .c[vl_covBed(.c, unique(peaks_rep[rep=="rep1"]))>0]
    .c <- .c[vl_covBed(.c, unique(peaks_rep[rep=="rep2"]))>0]
    .c$idx <- vl_collapseBed(.c, 
                             mingap = dist_cutoff, 
                             return_idx_only = T)
    .c[, c("start", "end"):= .(min(start), max(end)), idx]
    .c$idx <- NULL
    .c <- .c[, .SD[which.max(qValue)], .(seqnames, start, end)]
    fwrite(.c,
           filtered_peaks,
           sep= "\t",
           quote= F,
           na= ".",
           col.names= F)
  }
  print("DONE")
}, .(filtered_peaks, peaks_merge, dist_cutoff)]

#--------------------------------------------------------------#
# Merged_peaks
#--------------------------------------------------------------#
meta[!is.na(input), merged_file:= paste0("db/peaks/cutnrun_merged_peaks/", ChIP, "_merged_peaks", 
                                         ifelse(broad, ".broadPeak", ".narrowPeak"))]
meta[!is.na(input), {
  if(!file.exists(merged_file))
  {
    .c <- vl_importBed(unique(filtered_peaks))
    .c$idx <- vl_collapseBed(.c, 
                             mingap = dist_cutoff, 
                             return_idx_only = T)
    .c[, c("start", "end"):= .(min(start), max(end)), idx]
    .c$idx <- NULL
    .c <- .c[, .SD[which.max(qValue)], .(seqnames, start, end)]
    fwrite(.c,
           merged_file,
           sep= "\t",
           quote= F,
           na= ".",
           col.names= F)
  }
  print("DONE")
}, .(merged_file, dist_cutoff)]

#--------------------------------------------------------------#
# bw files
#--------------------------------------------------------------#
meta[!is.na(peaks_rep), bw_file:= paste0("db/bw/cutnrun/", ChIP, "_", cdition, "_", rep, ".bw")]
meta[!is.na(peaks_rep), {
  if(!file.exists(bw_file))
  {
    .g <- fread(gsub(ifelse(broad, "peaks.broadPeak$", "peaks.narrowPeak$"),
                     "treat_pileup.bdg" ,
                     peaks_rep),
                col.names = c("seqnames", "start", "end", "score"))
    .g[, start:= start+1]
    # Format GRanges
    .g <- GRanges(.g)
    BS <- BSgenome::getBSgenome("dm6")
    GenomeInfoDb::seqlevels(.g, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(BS)
    GenomeInfoDb::seqlengths(.g) <- GenomeInfoDb::seqlengths(BS)
    # save
    rtracklayer::export.bw(.g,
                           con= bw_file)
  }
  print("done")
}, .(peaks_rep, bw_file, broad)]
meta[!is.na(peaks_rep), bw_merge:= paste0("db/bw/cutnrun/", ChIP, "_", cdition, "_merge.bw")]
meta[!is.na(peaks_rep), {
  if(!file.exists(bw_merge))
  {
    vl_bw_merge(bw_file, 
                "dm6", 
                bins_width = 25L, 
                output = bw_merge)
  }
  print("done")
}, .(peaks_rep, bw_file, broad)]

#--------------------------------------------------------------#
# SAVE
#--------------------------------------------------------------#
fwrite(processed,
       "Rdata/processed_metadata_CUTNRUN.txt",
       na= NA)
