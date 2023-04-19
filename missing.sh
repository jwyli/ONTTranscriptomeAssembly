#!/bin/bash

refAnnotation=$1

grep -oE 'transcript_id "\S*"' $refAnnotation | cut -f2 -d'"'| sort | uniq > refTranscripts.tmp
grep -oE 'gene_id "\S*"' $refAnnotation | cut -f2 -d'"'| sort | uniq > refGenes.tmp

for assembler in flair stringtie.guided stringtie.free rb rattle
do
    grep -v 'novel\|associated_transcript' *$assembler*classification.txt| cut -f8 | sort | uniq > $assembler.transcripts.tmp
    grep -Fxv -f $assembler.transcripts.tmp refTranscripts.tmp > $assembler.transcripts.diff

    cut -f7 *$assembler*classification.txt | grep -v associated_gene | sort | uniq > $assembler.genes.tmp
    grep -Fxv -f $assembler.genes.tmp refGenes.tmp > $assembler.genes.diff
done
