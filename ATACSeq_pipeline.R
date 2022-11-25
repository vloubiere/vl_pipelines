setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(GenomicRanges)
require(readxl)

meta <- as.data.table(read_xlsx("Rdata/metadata_ATAC.xlsx"))
meta[, fq_file:= list.files("/mnt/f/_R_data/projects/epigenetic_cancer/db/fastq/ATAC_JJ_2018/", fq_basename, full.names = T), fq_basename]
meta[, bam_file:= paste0("/mnt/f/_R_data/projects/epigenetic_cancer/db/bam/ATAC/", gsub(".fastq.gz", ".bam", fq_basename)), fq_basename]
meta[, peaks_file:= paste0("db/peaks/ATAC/", gsub(".fastq.gz$", "", fq_basename), "_peaks.narrowPeak"), fq_basename]
meta[, bw_file:= paste0("db/bw/ATAC/", gsub(".fastq.gz", ".bw", fq_basename)), fq_basename]

meta[, {
  if(!file.exists(bam_file))
  {
    # bowtie 2
    sam_file <- gsub(".bam$", ".sam", bam_file)
    cmd <- paste0("bowtie2 -p 10 -x /mnt/d/_R_data/genomes/dm6/Sequence/Bowtie2Index/genome --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700")
    cmd <- paste0(cmd, " -U ", fq_file)
    cmd <- paste0(cmd, " -S ", sam_file)
    system(cmd)
    # sam to bam
    cmd <- paste0("/usr/bin/samtools view -@ 9 -b -q 30 ", sam_file, " -o ", bam_file)
    system(cmd)
    file.remove(sam_file)
  }
  # Peak calling
  if(!file.exists(peaks_file))
  {
    cmd <- paste0("/home/vloubiere/.local/bin/macs2 callpeak --keep-dup 1 -g dm --keep-dup 1 -f BAM -B --SPMR --nomodel --extsize 200",
                  " --outdir db/peaks/ATAC/ ",
                  " -t ", bam_file,
                  " -n ", gsub("_peaks.narrowPeak", "", basename(peaks_file)))
    system(cmd)
  }
  # Generate bw
  if(!file.exists(bw_file))
  {
    .g <- fread(gsub("peaks.narrowPeak$", "treat_pileup.bdg", peaks_file), col.names = c("seqnames", "start", "end", "score"))
    .g[, start:= start+1]
    # Format GRanges
    .g <- GRanges(.g)
    BS <- BSgenome::getBSgenome("dm6")
    GenomeInfoDb::seqlevels(.g) <- GenomeInfoDb::seqlevels(BS)
    GenomeInfoDb::seqlengths(.g) <- GenomeInfoDb::seqlengths(BS)
    # save
    rtracklayer::export.bw(.g,
                           con= bw_file)
  }
  print("DONE")
}, .(bam_file, peaks_file, bw_file)]

if(!file.exists("db/peaks/ATAC/ATAC_merged_peaks.narrowPeak"))
{
  cmd <- paste0("/home/vloubiere/.local/bin/macs2 callpeak --keep-dup 1 -g dm --keep-dup 1 -f BAM --nomodel --extsize 200",
                " --outdir db/peaks/ATAC/ ",
                " -t ", paste(meta$bam_file, collapse= " "),
                " -n ATAC_merged")
  system(cmd)
}

if(!file.exists("db/peaks/ATAC/ATAC_confident_peaks.narrowPeak"))
{
  conf <- vl_importBed("db/peaks/ATAC/ATAC_merged_peaks.narrowPeak")[signalValue>2 & qValue>2]
  rep_peaks <- rbindlist(lapply(meta$peaks_file, vl_importBed, extraCols= "narrowPeak"), idcol = T)
  conf <- conf[rep_peaks[conf, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N>=10]
  fwrite(conf[order(seqnames, start)],
         "db/peaks/ATAC/ATAC_confident_peaks.narrowPeak",
         col.names = F,
         sep= "\t",
         quote= F,
         na= NA)
}

if(!file.exists("db/bw/ATAC/ATAC_merged.bw"))
  vl_bw_merge(meta$bw_file, "dm6", bins_width = 25L, output = "db/bw/ATAC/ATAC_merged.bw")
