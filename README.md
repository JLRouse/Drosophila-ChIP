# Example of Drosophila-ChIP
Pipeline of ChIP-seq analysis in Drosophila utilising commonly used ChIP software

# Fastqc
To check sequencing data run through fastqc which should give an idea of the quality of the sequecning data

# Trimmomatic
To trim sequences that are reducing the quaility of the sequencing data, fastqc should again be run on the samples that have has some sequences trimmed.

Specific sequences can be trimmed using trimmomatic be creating FASTA files that can then be passed to the software. To create a FASTA file with a sequence to be trimmed open a word file via nano, then make sure the file starts with ">" and is described by .fa when saved. For example

```
>2R_input_repeated
TCTCTTCTCTTCTCTTCTCTTCTCTTCTCTTCTCTTCTCTTCTCTTCTC
>2R_input_repeated2
CAAGACAAGACAAGACAAGACAAGACAAGACAAGACAAGACAAGACAAGA
>2R_input_repeated3
AGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAG
```

To Trim using trimmomatic is it important to use both forward and reverse samples otherwise downstream mapping steps will not work. If this doesn't occur trimmomatic will trim bad quality sequences from one a single ended sequence and not its paired en alignment. This will lead to mapping software not able to pair the two sequences up
 
Example
```
java -jar ./trimmomatic-0.38-1/share/trimmomatic-0.38-1/trimmomatic.jar PE \
3_2S_4_1.fastq 3_2S_4_2.fastq 3_2S_4_1_trim.fastq 3_2S_4_1_unpaired.fastq \
3_2S_4_2_trim.fastq 3_2S_4_2_unpaired.fastq \
ILLUMINACLIP:./trimmomatic-0.38-1/share/trimmomatic-0.38-1/adapters/TruSeq3-SE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

Or to use own adaptors shown above
```
java -jar ./trimmomatic-0.38-1/share/trimmomatic-0.38-1/trimmomatic.jar PE \
2R_input_1.fastq 2R_input_2.fastq 2R_input_1_trim.fastq 2R_input_1_unpaired.fastq \
2R_input_2_trim.fastq 2R_input_2_unpaired.fastq \
ILLUMINACLIP:./trimmomatic-0.38-1/share/trimmomatic-0.38-1/adapters/2R_adapter.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

# Align using bowtie 2
Can align using bowtie 2. First, bowtie 2 requires a reference genome built from fasta file. In this case I am using the Drosophila melanogaster genome and saving is as "bt2reference"

```
bowtie2-build /nobackup/fbsjrou/reference_annotation/Drosophila_melanogaster.BDGP6.dna.chromosome.2L.fa,\
/nobackup/fbsjrou/reference_annotation/Drosophila_melanogaster.BDGP6.dna.chromosome.2R.fa,\
/nobackup/fbsjrou/reference_annotation/Drosophila_melanogaster.BDGP6.dna.chromosome.3L.fa,\
/nobackup/fbsjrou/reference_annotation/Drosophila_melanogaster.BDGP6.dna.chromosome.3R.fa,\
/nobackup/fbsjrou/reference_annotation/Drosophila_melanogaster.BDGP6.dna.chromosome.4.fa,\
/nobackup/fbsjrou/reference_annotation/Drosophila_melanogaster.BDGP6.dna.chromosome.X.fa,\
/nobackup/fbsjrou/reference_annotation/Drosophila_melanogaster.BDGP6.dna.chromosome.Y.fa bt2reference
```

Should then align own files (which are trimmed) to this reference genome. This is easiest if the files are all in the same directory. For example, here both the forward (1) and reverse (2) reads for 3_2R_4 and 3_2S_4 are aligned. The defined .sam files are output with the names 3_2R_4 or 3_2S_4

```
bowtie2 - p 8  -x bt2reference -1 /nobackup/fbsjrou/ChIP_gz/trimmed/3_2R_4_1_trim.fastq -2 /nobackup/fbsjrou/ChIP_gz/trimmed/3_2R_4_2_trim.fastq -S 3_2R_4.sam
bowtie2 -p 8 -x bt2reference -1 /nobackup/fbsjrou/ChIP_gz/trimmed/3_2S_4_1_trim.fastq -2 /nobackup/fbsjrou/ChIP_gz/trimmed/3_2S_4_2_trim.fastq -S 3_2S_4.sam
```

Aligning should take place for all samples of interest, including input samples. The next step will start to seperate output ChIpped samples and input samples

# Call peaks using MACS2
Once ChIP sequences are aligned, need to call peaks using Macs2. This program takes input samples and compares the number of reads to ChIped samples. Where ChIped samples have a much larger number of reads the software will call a peak. It is necessary to change the parameters of Macs2 to define the best peak calling parameters, this is especailly important when investgiating histone marks.

For example
```
macs2 callpeak -t /nobackup/fbsjrou/ChIP_gz/mapping/3_2R_4_sorted.bam \
-c /nobackup/fbsjrou/ChIP_gz/mapping/2R_input.bam \
-B --nomodel --extsize 147 --SPMR --keep-dup auto -f BAMPE -g dm -n 3_2R_4
```

To detail what is written above

callpeak from the experimental samples, with full pathway to .bam files
```
macs2 callpeak -t /nobackup/fbsjrou/ChIP_gz/mapping/3_2R_4_sorted.bam \
```
Specify control samples with full pathway to .bam files, make sure these pair up with the sample above
```
-c /nobackup/fbsjrou/ChIP_gz/mapping/2R_input.bam \
```
Stores the fragment pileup, control lambda and bedGraph files, important for visualisation
```
-B
```
Specify no model to be shown, which reduces the disk space needed
```
--nomodel
```
Specifies the size that the genome is cut into, use 147 in this case as that is how many bases are wrapped around a histone
```
--extsize 147
```
Generates bdg files
```
--SPMR
```
Keep duplicates of peaks
```
--keep-dup auto 
```
Specify format, this is .bam files paired end sequencing 
```
-f BAMPE
```
Specifies genome size, dm stands for Drosophila melanogaster
```
-g dm 
```
Specifies what you name all this. All files will get this name just with different file extensions. 
```
-n 3_2S_4
```

# Calling differential peaks using ChIPComp
To call differential peaks with ChIPComp you will need a peak file (.bed, provided by Macs2) and control file (.bam, provided by bowtie 2) and a sample file (.bam, provided by bowtie 2). All file locations should be defined in a table that can then be passed to the ChIPComp software. I find this is most easily done with a csv file
For example ("ChIPsetup4.csv" below)

SampleID	condition	factor	ipReads	ctReads	peaks
1	single	H3K4me3	3_2S_4_sorted.bam	2S_input_sorted.bam	3_2S_4_peaks.bed
2	single	H3K4me3	4_2S_4_sorted.bam	2S_input_sorted.bam	4_2S_4_peaks.bed
3	single	H3K4me3	5_2S_4_sorted.bam	2S_input_sorted.bam	5_2S_4_peaks.bed
4	rival	H3K4me3	3_2R_4_sorted.bam	2R_input_sorted.bam	3_2R_4_peaks.bed
5	rival	H3K4me3	4_2R_4_sorted.bam	2R_input_sorted.bam	4_2R_4_peaks.bed
6	rival	H3K4me3	5_2R_4_sorted.bam	2R_input_sorted.bam	5_2R_4_peaks.bed

```R
library("ChIPComp")
confs=makeConf("ChIPsetup4.csv")
conf=confs$conf
design=confs$design
countSet=makeCountSet(conf,design,filetype="bam",binsize=50, species="other", mva.span=c(1000))
pdf("2_hours_4_1000.pdf")
plot(countSet)
dev.off()
countSet=ChIPComp(countSet)
print(countSet)
write.csv(countSet$db, "2_hours_4_1000.csv")
```

# Defining peaks with ChIPseeker
ChIPseeker allows the annotation of peaks and some really cool graphical visulisations that come in handy when organising your ChIP-seq data
Very simply, to annotate diffeential peaks, ChIPseeker will compare genomic regions to gff3 or gtf files and assign an annotation (gene and feature) to a differential peak.

```r
library("rlang")
library("ChIPseeker")
library("GenomicFeatures")
library("rtracklayer")
library("RSQLite")
txdb <- makeTxDbFromGFF("Drosophila_melanogaster.BDGP6.22.96.chr.gff3", format="gff3")
promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
DE = "H3K4 diffregulated.txt"
peak <- readPeakFile(DE, as="GRanges", header=TRUE)
tagMatrixDE <- getTagMatrix(peak, windows=promoter)
pdf("diffregpromoter.pdf")
tagHeatmap(tagMatrixDE, xlim=c(-2000, 2000), color="red")
dev.off()
peakAnnoListDE <- annotatePeak(DE, tssRegion=c(-2000, 2000), TxDb=txdb, addFlankGeneInfo = TRUE, flankDistance = 5000, overlap="all")
pdf("Allplotdiffreg.pdf")
plotAnnoPie(peakAnnoListDE)
vennpie(peakAnnoListDE)
upsetplot(peakAnnoListDE)
upsetplot(peakAnnoListDE, vennpie=TRUE)
dev.off()
write.csv(peakAnnoListDE, "H3K4diffreg.csv")
```

# Clustering peaks
Off the back of differential peak calling with ChIPComp peaks can be clustered to provide an overview of what is occuring in the genome. It is exploratory but can be very powerful. This should be thought of in a similar way to clustering analysis in RNA-seq using normalised counts but instead using normalised peak regions.


