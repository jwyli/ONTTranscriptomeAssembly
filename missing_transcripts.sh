#!/bin/bash

refAnnotation=$1

grep transcript $refAnnotation | cut -f9 | cut -f4 -d'"'|uniq >> refGenes.tmp

for assembler in flair stringtie.guided stringtie.free rb rattle
do
    grep -v 'novel\|associated_transcript' *$assembler*classification.txt| cut -f8 | sort | uniq > $assembler.tmp
    grep -Fxv -f $assembler.tmp refGenes.tmp > $assembler.diff
done

rm *.tmp
