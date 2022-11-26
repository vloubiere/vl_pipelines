setwd("/groups/stark/vloubiere/projects/insulator_tomas/")
require(vlfunctions)
require(readxl)
require(data.table)
require(parallel)
require(Rsubread)
require(seqinr)
require(rtracklayer)

# Metadata (select relevant tracks)
CTCF <- fread("db/metadata/SraRunTable_CTCF_bingRen.txt")[Genotype=="Wild-type" & `Assay Type`=="ChIP-Seq"]
Ab <- data.table(chip_antibody=c("ab992", "ab70303", "Active Motif\\,39685", "Millipore\\, 04-745", "none"),
                 cdition= "WT",
                 ChIP= c("Rad21", "CTCF", "H3K27Ac", "H3K4me3", "Input"))
CTCF <- CTCF[Ab, on= "chip_antibody"]
MAZ <- fread("db/metadata/SraRunTable_MAZ_Reinberg.txt")[Genotype=="wild type" & `Assay Type`=="ChIP-Seq" & source_name=="Embryonic stem cells"]
Ab <- data.table(Run= c("SRR17298884", "SRR17298890", "SRR17298874"),
                 cdition= "WT",
                 ChIP= c("MAZ", "MAZ", "Input"))
MAZ <- MAZ[Ab, on= "Run"]
meta <- rbindlist(list(Ren= CTCF, Reinberg= MAZ), idcol = "lab", fill= T)
meta <- meta[, .(Run, LibraryLayout, lab, cdition, ChIP)]
meta[, rep:= paste0("rep", seq(.N)), .(lab, ChIP)]
meta[, sample:= paste0(lab, "_", cdition, "_", ChIP, "_", rep)]
# Set working directory
meta[, cmd:= paste(c(paste("cd", getwd()),
                     "module load build-env/2020", 
                     "sra-toolkit/2.9.6-1-centos_linux64",
                     "module load trim_galore/0.6.2-foss-2018b-python-3.6.6",
                     "module load bowtie2/2.3.5.1-foss-2018b", 
                     "module load samtools/1.9-foss-2018b",
                     "module load macs2/2.2.5-foss-2018b-python-3.6.6"), 
                   collapse= "; ")]

#--------------------------------------------------------------#
# Download fastq files using SRR
#--------------------------------------------------------------#
meta[, dl_cmd:= {
  if(length(list.files("db/fastq/", paste0(Run, ".*.fastq.gz")))==0)
    paste0("fastq-dump -A ", 
           Run, 
           " --split-files --gzip -O /groups/stark/vloubiere/projects/insulator_tomas/db/fastq/")
}, Run]

#--------------------------------------------------------------#
# Trim reads
#--------------------------------------------------------------#
meta[, fq1:= paste0("db/fastq/", Run, "_1.fastq.gz"), Run]
meta[, fq2:= paste0("db/fastq/", Run, "_2.fastq.gz"), Run]
meta[LibraryLayout=="PAIRED", fq1_trim:= gsub(".fastq.gz$", "_val_1.fq.gz", fq1)]
meta[LibraryLayout=="PAIRED", fq2_trim:= gsub(".fastq.gz$", "_val_2.fq.gz", fq2)]
meta[LibraryLayout=="SINGLE", fq1_trim:= gsub(".fastq.gz$", "_trimmed.fq.gz", fq1)]
meta[, trim_cmd:= {
  outdir <- dirname(fq1)
  if(!file.exists(fq1_trim) && LibraryLayout=="SINGLE")
  {
    paste0("trim_galore --gzip -o ", outdir, "/ ", fq1)
  }else if((!file.exists(fq1_trim) | !file.exists(fq2_trim)) && LibraryLayout=="PAIRED")
    paste0("trim_galore --paired --gzip -o ", outdir, "/ ", fq1, " ", fq2)
}, .(fq1, fq2, fq1_trim, fq2_trim, LibraryLayout)]

#--------------------------------------------------------------#
# Alignment
#--------------------------------------------------------------#
meta[, bam:= paste0("db/bam/", lab, "_", cdition, "_", ChIP, "_", rep, ".bam")]
meta[, align_cmd:= {
  if(!file.exists(bam))
  {
    # bowtie 2
    sam_file <- gsub(".bam$", ".sam", bam)
    aln <- "bowtie2 -p 10 -x /groups/stark/vloubiere/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700"
    if(LibraryLayout=="SINGLE")
      aln <- paste0(aln, " -U ", fq1_trim) else if(LibraryLayout=="PAIRED")
        aln <- paste0(aln, " -1 ", fq1_trim, " -2 ", fq2_trim)
    aln <- paste0(aln, " -S ", sam_file)
    # sam to bam
    filter <- paste0("samtools view -@ 9 -b -q 30 ", sam_file, " -o ", bam)
    # Remove sam file
    clean <- paste0("rm ", sam_file)
    cmd <- c(aln, filter, clean)
  }
  paste(cmd, collapse= "; ")
}, .(bam, LibraryLayout)]

#--------------------------------------------------------------#
# Peak calling separated replicates
#--------------------------------------------------------------#
meta[lab=="Ren" & ChIP!="Input", input_bam:= meta[lab=="Ren" & ChIP=="Input"][.BY, bam, on= "rep"], rep]
meta[lab=="Ren" & ChIP!="Input" & is.na(input_bam), input_bam:= meta[lab=="Ren" & ChIP=="Input" & rep=="rep2", bam]]
meta[lab=="Reinberg" & ChIP!="Input", input_bam:= meta[lab=="Reinberg" & ChIP=="Input" & rep=="rep1", bam]]
meta[!is.na(input_bam), peaks_rep:= paste0("db/peaks/", sample)]
meta[!is.na(input_bam), peakRep_cmd:= {
  if(!file.exists(paste0(peaks_rep, "_peaks.narrowPeak")))
  {
    paste0("macs2 callpeak -g mm --keep-dup 1 --SPMR --outdir ", paste0(dirname(peaks_rep), "/"),
           " --nomodel --extsize 200 ",
           " -t ", paste(bam, collapse= " "), 
           " -c ", paste(input_bam, collapse= " "),
           " -n ", basename(peaks_rep))
  }
}, peaks_rep]

#--------------------------------------------------------------#
# Peak calling merged replicates
#--------------------------------------------------------------#
meta[!is.na(input_bam), peaks_merge:= gsub("^(.*)_rep.*", "\\1_merge", peaks_rep)]
meta[!is.na(input_bam), peakMerge_cmd:= {
  if(!file.exists(paste0(peaks_merge, "_peaks.narrowPeak")))
  {
    paste0("macs2 callpeak -g mm --keep-dup 1 -B --SPMR --outdir ", paste0(dirname(peaks_merge), "/"),
           " --nomodel --extsize 200 ",
           " -t ", paste(bam, collapse= " "),
           " -c ", paste(input_bam, collapse= " "),
           " -n ", basename(peaks_merge))
  }
}, peaks_merge]
meta[, final:= paste(na.omit(unique(unlist(.SD))), collapse = "; "), .(lab, cdition, ChIP), .SDcols= patterns("cmd$")]
# meta[, vl_bsub(final, o = "log/", cores = 8, m = 20, t= '08:00:00'), final]

#--------------------------------------------------------------#
# bw files
#--------------------------------------------------------------#
meta[, bw_file:= paste0("db/bw/", lab, "_", cdition, "_", ChIP, "_merge.bw")]
meta[!is.na(peaks_merge), {
  if(!file.exists(bw_file))
  {
    .g <- fread(paste0(peaks_merge, "_treat_pileup.bdg"),
                col.names = c("seqnames", "start", "end", "score"))
    .g[, start:= start+1]
    # Format GRanges
    .g <- GRanges(.g)
    BS <- BSgenome::getBSgenome("mm10")
    GenomeInfoDb::seqlevels(.g, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(BS)
    GenomeInfoDb::seqlengths(.g) <- GenomeInfoDb::seqlengths(BS)
    # save
    rtracklayer::export.bw(.g,
                           con= bw_file)
  }
  print("done")
}, .(peaks_merge, bw_file)]

#--------------------------------------------------------------#
# Confident peaks
#--------------------------------------------------------------#
meta[!is.na(peaks_merge), confident_peaks:= paste0(peaks_merge, "_confident.narrowPeak"), peaks_merge]
meta[!is.na(peaks_merge), {
  if(!file.exists(confident_peaks))
  {
    .c <- vl_importBed(paste0(peaks_merge, "_peaks.narrowPeak"))
    ov <- lapply(peaks_rep, function(x) vl_covBed(.c, paste0(x, "_peaks.narrowPeak")))
    .c <- .c[apply(as.data.table(ov), 1, function(x) all(x>0))]
    rtracklayer::export(.c, 
                        confident_peaks)  
  }
  .SD
}, .(peaks_merge, confident_peaks)]

#--------------------------------------------------------------#
# Screenshot
#--------------------------------------------------------------#
tracks <- c("db/bw/Reinberg_WT_MAZ_merge.bw",
            "db/peaks/Reinberg_WT_MAZ_merge_confident.narrowPeak",
            "db/bw/Ren_WT_CTCF_merge.bw",
            "db/peaks/Ren_WT_CTCF_merge_confident.narrowPeak",
            "db/bw/Ren_WT_H3K27Ac_merge.bw",
            "db/peaks/Ren_WT_H3K27Ac_merge_confident.narrowPeak",
            "db/bw/Ren_WT_H3K4me3_merge.bw",
            "db/peaks/Ren_WT_H3K4me3_merge_confident.narrowPeak",
            "db/bw/Ren_WT_Rad21_merge.bw",
            "db/peaks/Ren_WT_Rad21_merge_confident.narrowPeak")
pdf("pdf/screenshot_ChIPSeq.pdf", height = 5)
vl_screenshot(data.table(seqnames= "chr3", start= 50e6, end= 53e6),
              tracks = tracks,
              names= make.unique(unlist(tstrsplit(tracks, "_", keep= 3))))
dev.off()
