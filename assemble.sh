#!/bin/bash

homedir="/hdd1/home/p21_wyli/6010/"
baseForName=$1 
readsFq=$2 #full path
refGenome=$3 #full path
refAnnotation=$4 #full path
# dataPath=$(dirname "${readsFq}") ; dataFile=$(basename "${readsFq}")

exec > >(tee "${baseForName}.assembly.out") 2>&1 # save all subsequent output
set -x #echo on

# perform transcriptome assembly
cd ${homedir}/02.assembly
mkdir $baseForName
cd $baseForName 

# run flair
#check runtimes
echo "starting flair analysis for ${baseForName}"
starttime=`date +%s`

# run flair align
# usage: flair align -g genome.fa -r <reads.fq>|<reads.fa> [options]
flair align -g $refGenome -r $readsFq --output "${baseForName}.flair.aligned" --nvrna -t 48

# run flair correct
# usage: flair correct -q query.bed12 [-f annotation.gtf]|[-j introns.tab] -g genome.fa [options]
flair correct -q "${baseForName}.flair.aligned.bed" -f $refAnnotation -g $refGenome \
--output "${baseForName}.flair" --nvrna -t 48

# run flair collapse
# usage: flair collapse -g genome.fa -q <query.bed> -r <reads.fq>/<reads.fa> [options]
flair collapse -g $refGenome -q "${baseForName}.flair_all_corrected.bed" -r $readsFq --gtf $refAnnotation \
-t 48 --output "${baseForName}.flair.collapse"

#check runtimes
endtime=`date +%s`
runtime=$((endtime-starttime))
echo "finished flair analysis for ${baseForName} in ${runtime} seconds"

#run stringtie
#check runtimes
echo "starting stringtie analysis for ${baseForName}"
starttime=`date +%s`

#Input: spliced RNA-seq read alignments (SAM, BAM or CRAM file sorted by coordinate). 
#The -L option must be used when the input alignment file contains (sorted) spliced alignments of long read RNA-seq or cDNA reads. 
#Such alignments can be produced by `minimap2` with the `-ax splice` option, which also generates the necessary `ts` tag to indicate the transcription strand. 
#As mentioned above such `minimap2` alignment files must be first position-sorted by before they can be processed by StringTie. 
#mm2 github manual suggest: minimap2 -ax splice -uf -k14 ref.fa direct-rna.fq > aln.sam
#? --splice-flank=no if sirv? remember not to add for other samples! 
# https://github.com/lh3/minimap2/blob/master/cookbook.md#map-direct-rna

#-a 	Generate CIGAR and output alignments in the SAM format
#-x = preset values; put in front so that the subsequent values can override settings if needed
minimap2 -t 48 -ax splice -uf -k14 $refGenome $readsFq \
| samtools sort -T ./tmp -O bam -o "${baseForName}.stringtie.aln.bam"

#sgnex says:
# use “-ax splice --junc-bed” for genomic alignments using the junction bed file to correct splicing junctions, 
# use “-ax map-ont” for transcriptomic alignments. 
# For direct RNA-Seq runs the additional parameters “-k14” and “-uf” are used as recommended

# https://github.com/lh3/minimap2#map-long-splice says:
# For Iso-seq, Direct RNA-seq and tranditional full-length cDNAs, 
# it would be desired to apply -u f to force minimap2 to consider the forward transcript strand only. 
# This speeds up alignment with slight improvement to accuracy. 
# For noisy Nanopore Direct RNA-seq reads, it is recommended to use a smaller k-mer size for increased sensitivity to the first or the last exons.

samtools index "${baseForName}.stringtie.aln.bam"

#A reference annotation file in GTF or GFF3 format can be provided to StringTie using the -G option 
#used as 'guides' for the assembly process and help improve the transcript structure recovery for those transcripts.
#rec if model organism eg. human

#usage: stringtie [-o <output.gtf>] [other_options] <read_alignments.bam>
#without guide annotation
# -c: Sets the minimum read coverage allowed for the predicted transcripts. 
# A transcript with a lower coverage than this value is not shown in the output. Default: 1 
# -s: minimum read coverage for single-exon transcripts. default = 4.75 
#not sure why rnabloom supp set as 3??
stringtie -p 48 -L -c 3 -s 3 \
-o "${baseForName}.stringtie.free.assembly.gtf" "${baseForName}.stringtie.aln.bam"

#with guide annotation
stringtie -p 48 -L -c 3 -s 3 -G $refAnnotation \
-o "${baseForName}.stringtie.guided.assembly.gtf" "${baseForName}.stringtie.aln.bam"

#use gffread to generate a FASTA file with the DNA sequences for all transcripts in a GFF file
#usage: gffread -w transcripts.fa -g /path/to/genome.fa transcripts.gtf
gffread -w "${baseForName}.stringtie.free.assembly.fa" -g $refGenome "${baseForName}.stringtie.free.assembly.gtf"
gffread -w "${baseForName}.stringtie.guided.assembly.fa" -g $refGenome "${baseForName}.stringtie.guided.assembly.gtf"

# check runtimes
endtime=`date +%s`
runtime=$((endtime-starttime))
echo "finished stringtie analysis for ${baseForName} in ${runtime} seconds for both with and without guide gtf."

# run rnabloom2
#check runtimes
echo "starting rnabloom2 analysis for ${baseForName}"
starttime=`date +%s`

#usage: [java -jar RNA-Bloom.jar]|rnabloom -long LONG.fastq -stranded -t THREADS -outdir OUTDIR
rnabloom -long $readsFq -stranded -t 48 -outdir . -n "${baseForName}.rb"

# check runtimes
endtime=`date +%s`
runtime=$((endtime-starttime))
echo "finished rnabloom2 analysis for ${baseForName} in ${runtime} seconds"


