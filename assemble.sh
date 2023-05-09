#!/bin/bash

homedir="/hdd1/home/p21_wyli/6010/"
baseForName=$1 
readsFq=$2 #full path
refGenome=$3 #full path
refAnnotation=$4 #full path

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
starttime_map=`date +%s`

#Input: spliced RNA-seq read alignments (SAM, BAM or CRAM file sorted by coordinate). 
#-a 	Generate CIGAR and output alignments in the SAM format -x = preset values; put in front so that the subsequent values can override settings if needed
minimap2 -t 48 -ax splice -uf -k14 $refGenome $readsFq \
| samtools sort -T ./tmp -O bam -o "${baseForName}.stringtie.aln.bam"

# For noisy Nanopore Direct RNA-seq reads, it is recommended to use a smaller k-mer size for increased sensitivity to the first or the last exons.
samtools index "${baseForName}.stringtie.aln.bam"

endtime_map=`date +%s`
runtime_map=$((endtime_map-starttime_map))

#without guided annotation
starttime=`date +%s`

stringtie -p 48 -L -o "${baseForName}.stringtie.free.assembly.gtf" "${baseForName}.stringtie.aln.bam"

#use gffread to generate a FASTA file with the DNA sequences for all transcripts in a GFF file
#usage: gffread -w transcripts.fa -g /path/to/genome.fa transcripts.gtf
gffread -w "${baseForName}.stringtie.free.assembly.fa" -g $refGenome "${baseForName}.stringtie.free.assembly.gtf"

endtime=`date +%s`
runtime=$((endtime-starttime+runtime_map))
echo "finished stringtie.free analysis for ${baseForName} in ${runtime} seconds."

#with guide annotation
starttime=`date +%s`

stringtie -p 48 -L -G $refAnnotation -o "${baseForName}.stringtie.guided.assembly.gtf" "${baseForName}.stringtie.aln.bam"
gffread -w "${baseForName}.stringtie.guided.assembly.fa" -g $refGenome "${baseForName}.stringtie.guided.assembly.gtf"

# check runtimes
endtime=`date +%s`
runtime=$((endtime-starttime+runtime_map))
echo "finished stringtie.guided analysis for ${baseForName} in ${runtime} seconds."

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
