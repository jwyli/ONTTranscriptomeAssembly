#!/bin/bash

refAnnotation=$1

grep -oE 'transcript_id "\S*"' $refAnnotation | cut -f2 -d'"'| sort | uniq > refGenes.tmp

for assembler in flair stringtie.guided stringtie.free rb rattle
do
    grep -v 'novel\|associated_transcript' *$assembler*classification.txt| cut -f8 | sort | uniq > $assembler.tmp
    grep -Fxv -f $assembler.tmp refGenes.tmp > $assembler.diff
done
